function plotFFT(signal, Fs)
% plotFFT - Plots the FFT of a signal
% Inputs:
%   signal - input time-domain signal (vector)
%   Fs     - sampling frequency in Hz
%
% Usage:
%   plotFFT(signal, Fs)

    % Ensure the signal is a column vector
    signal = signal(:);
    
    % Length of the signal
    N = length(signal);
    
    % Compute the FFT
    Y = fft(signal);
    
    % Compute the two-sided spectrum
    P2 = abs(Y / N);
    
    % Compute the single-sided spectrum
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Frequency vector
    f = Fs * (0:floor(N/2)) / N;
    
    % Plot
    figure;
    plot(f, P1);
    title('Single-Sided Amplitude Spectrum of Signal');
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
    grid on;
end
