function enhanced_signal = rls_filter(noisy_speech, external_noise, fs, num_taps, lambda, R, notch, tones)
    delta = 0.001;
    w_rls = zeros(num_taps, 1);
    P = (1/delta) * eye(num_taps);

    filtered_buffer = zeros(size(external_noise));
    enhanced_signal = zeros(size(noisy_speech));

    %construct notch filter if partial suppression is enabled
    if notch == 1
        [b_all, a_all] = design_notch_filter(tones, fs, R);
        y_state = zeros(length(tones), 2);
        x_state = zeros(length(tones), 2);
    end

    for n = 1:length(noisy_speech)
        % apply notch on reference is partial mode is on.
        if notch == 1
            [filtered_val, y_state, x_state] = apply_cascade_notch(external_noise(n), y_state, x_state, b_all, a_all);
            filtered_buffer(n) = filtered_val;
        else
            filtered_buffer(n) = external_noise(n); % no notch reqd. for full suppr.
        end
        
        % RLS equations
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
