function filtered_signal = ideal_notch_filter(input_signal, fs, notchFreq, notchWidth)
    N = length(input_signal);
    X = fft(input_signal);
    f_fft = (0:N-1) * (fs/N);

    notchBins = abs(f_fft - notchFreq) <= notchWidth/2;

    negNotchFreq = fs - notchFreq;
    negNotchBins = abs(f_fft - negNotchFreq) <= notchWidth/2;

    X(notchBins) = 0;
    X(negNotchBins) = 0; 

    filtered_signal = real(ifft(X));
end