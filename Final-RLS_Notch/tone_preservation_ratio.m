function TPI = tone_preservation_ratio(noisy_signal, filtered_signal, tonal_freq, fs)
    N = length(noisy_signal);
    f = (0:N-1)*(fs/N);
    
    fft_noisy = abs(fft(noisy_signal));
    fft_filtered = abs(fft(filtered_signal));

    %the closest frequency to the tone in "f"
    [~, idx] = min(abs(f - tonal_freq));
    TPI = fft_filtered(idx) / fft_noisy(idx);
end