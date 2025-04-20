function enhanced_signal = rls_filter(noisy_speech, external_noise, num_taps, lambda, notch)
    delta = 0.001; 
    w_rls = zeros(num_taps, 1); 
    P = (1/delta) * eye(num_taps); 
    fs = 44100;

    if (notch == 1)
        notchFreq = 1000;
        Q = 35;
        w0 = 2 * pi * notchFreq / fs;
        alpha = sin(w0) / (2 * Q);

        b0 = 1;
        b1 = -2 * cos(w0);
        b2 = 1;
        a0 = 1 + alpha;
        a1 = -2 * cos(w0);
        a2 = 1 - alpha;

        b = [b0, b1, b2] / a0;
        a = [1, a1 / a0, a2 / a0];
    end

    % Initialize filtered buffer
    filtered_buffer = zeros(size(external_noise));
    enhanced_signal = zeros(size(noisy_speech));

    for n = 1:length(noisy_speech)
        % Only compute filtered_buffer(n) once, based on prior samples
        if notch == 1
            x0 = external_noise(n);
            if n > 1
                x1 = external_noise(n-1);
                y1 = filtered_buffer(n-1);
            else
                x1 = 0;
                y1 = 0;
            end
            if n > 2
                x2 = external_noise(n-2);
                y2 = filtered_buffer(n-2);
            else
                x2 = 0;
                y2 = 0;
            end

            filtered_buffer(n) = b(1)*x0 + b(2)*x1 + b(3)*x2 - a(2)*y1 - a(3)*y2;
        else
            filtered_buffer(n) = external_noise(n); % no filtering
        end

        % Start RLS only after num_taps samples
        if n >= num_taps
            x_vec = filtered_buffer(n - num_taps + 1 : n);
            x_vec = flip(x_vec);

            d = noisy_speech(n);

            K = (P * x_vec) / (lambda + x_vec' * P * x_vec);
            v_hat = w_rls' * x_vec;
            e = d - v_hat;

            w_rls = w_rls + K * e;
            P = (P - K * x_vec' * P) / lambda;

            enhanced_signal(n) = e;
        end
    end
end
