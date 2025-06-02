function plot_spectrogram(audio_signal, Fs, type)
    window_size = 1024;  
    overlap = window_size / 2;  
    nfft = 2^nextpow2(window_size);  

    [S, F, T] = spectrogram(audio_signal, window_size, overlap, nfft, Fs);
    
    S_dB = 20 * log10(abs(S));  
    
    min_dB = -120;
    max_dB = 40;
    
    figure;
    imagesc(T, F, S_dB);
    axis xy;
    
    colormap jet;  
    colorbar;  
    clim([min_dB max_dB]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    if (type == 0)
        title('Spectrogram of Clean Signal (dB)');
    elseif (type == 1)
        title('Spectrogram of Full-Suppression Mode Output (dB)');
    elseif (type == 2)
        title('Spectrogram of Partial-Suppression Output (dB)');
    end
end
