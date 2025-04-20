% %function filtered_signal = practical_notch_filter(signal, fs, notchFreq, Q)
%     wo = notchFreq / (fs / 2); 
%     bw = wo / Q; 
%     [b, a] = iirnotch(wo, bw);
%     filtered_signal = filter(b, a, signal);
% end


function filtered_signal = practical_notch_filter(signal, fs, notchFreq, Q)
    % Normalize the notch frequency
    wo = 2 * pi * notchFreq / fs;

    % Compute the filter coefficients manually
    alpha = sin(wo) / (2 * Q);

    b0 =  1;
    b1 = -2 * cos(wo);
    b2 =  1;
    a0 =  1;
    a1 = -2 * cos(wo) * sqrt(1 - alpha);
    a2 =  1 - alpha;

    % Normalize coefficients
    b0 = b0 / a0;
    b1 = b1 / a0;
    b2 = b2 / a0;
    a1 = a1 / a0;
    a2 = a2 / a0;

    % Preallocate output
    filtered_signal = zeros(size(signal));

    % Initialize past samples
    x1 = 0; x2 = 0;
    y1 = 0; y2 = 0;

    % Apply the filtering using the difference equation
    for n = 1:length(signal)
        x0 = signal(n);
        y0 = b0*x0 + b1*x1 + b2*x2 - a1*y1 - a2*y2;

        filtered_signal(n) = y0;

        % Update past samples
        x2 = x1;
        x1 = x0;
        y2 = y1;
        y1 = y0;
    end
end
