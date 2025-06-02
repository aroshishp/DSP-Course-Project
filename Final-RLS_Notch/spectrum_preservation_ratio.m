function ratio = spectrum_preservation_ratio(clean_signal, filtered_signal, tonal_freq, fs)
    N = length(clean_signal);
    f = (0:N-1)*(fs/N);
    
    % limit according to Nyquist frequency
    half_N = ceil(N/2);
    f_half = f(1:half_N);
    
    % non-tonal regions
    mask = (f_half < tonal_freq - 20) | (f_half > tonal_freq + 20);

    fft_clean = abs(fft(clean_signal));
    fft_filtered = abs(fft(filtered_signal));
    
    fft_clean_masked = fft_clean(1:half_N);
    fft_filtered_masked = fft_filtered(1:half_N);
    
    ratio = sum(fft_filtered_masked(mask).^2) / sum(fft_clean_masked(mask).^2);
end