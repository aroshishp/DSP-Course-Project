function [notched, removed] = split_notch_filter(signal, fs, notchFreq, Q)
    wo = 2 * pi * notchFreq / fs;
    alpha = sin(wo) / (2 * Q);

    b0 = 1;
    b1 = -2 * cos(wo);
    b2 = 1;
    a0 = 1 + alpha;
    a1 = -2 * cos(wo);
    a2 = 1 - alpha;

    b = [b0, b1, b2] / a0;
    a = [1, a1 / a0, a2 / a0];

    notched = zeros(size(signal));
    removed = zeros(size(signal));
    x1 = 0; x2 = 0; y1 = 0; y2 = 0; yn1 = 0; yn2 = 0;

    for n = 1:length(signal)
        x0 = signal(n);
        y0 = b(1)*x0 + b(2)*x1 + b(3)*x2 - a(2)*y1 - a(3)*y2;
        rn = x0 - y0;

        notched(n) = y0;
        removed(n) = rn;

        x2 = x1; x1 = x0;
        y2 = y1; y1 = y0;
    end
end
