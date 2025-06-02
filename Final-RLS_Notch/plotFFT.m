function plotFFT(signal, Fs, type)
    signal = signal(:);
    N = length(signal);

    Y = fft(signal);
    P2 = abs(Y / N);
    
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs * (0:floor(N/2)) / N;
    
    figure;
    plot(f, P1);
    if (type == 1)
        title('Spectrum of Partial Suppression Output');
    else
        title('Spectrum of Signal');
    end
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
    grid on;
end
