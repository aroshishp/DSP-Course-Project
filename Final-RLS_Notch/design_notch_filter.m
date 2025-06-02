function [b_all, a_all] = design_notch_filter(tones, fs, R)
    num_tones = length(tones);
    b_all = zeros(num_tones, 3);
    a_all = zeros(num_tones, 3);

    for k = 1:num_tones
        w0 = 2 * pi * tones(k) / fs;

        b = [1, -2 * cos(w0), 1];
        a = [1, -2 * R * cos(w0), R^2];

        b_all(k,:) = b;
        a_all(k,:) = a;
    end
end
