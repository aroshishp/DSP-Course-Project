% Load signals
clean = load('clean_speech.txt');
noise = load('external_noise.txt');
noisy = load('noisy_speech.txt');

% Standardize signals
clean = clean(:); noise = noise(:); noisy = noisy(:);
min_len = min([length(clean), length(noise), length(noisy)]);
clean = clean(1:min_len);
noise = noise(1:min_len);
noisy = noisy(1:min_len);

fs = 44100;  % Sample rate
L = min_len;
t = (0:L-1)/fs;
f = fs*(0:(L/2))/L;

%% Initial SNR Calculation
original_noise = noisy - clean;
snr_before = 10*log10(sum(clean.^2)/sum(original_noise.^2));

%% ====== ADAPTIVE NOTCH FILTERING ======
% Frequency analysis
Y_noisy = fft(noisy);
spectrum = abs(Y_noisy(1:floor(L/2)+1));

% Find tonal components
[pks,locs] = findpeaks(spectrum,...
    'MinPeakHeight', 0.2*max(spectrum),...
    'MinPeakDistance', 50);

% Apply notch filters
bw = 25;
notch_filtered = noisy;
for k = 1:length(locs)
    freq = f(locs(k));
    if freq < 0.99*(fs/2)
        wo = freq/(fs/2);
        bw_norm = bw/(fs/2);
        [b,a] = iirnotch(wo, bw_norm);
        notch_filtered = filter(b,a,notch_filtered);
    end
end

% SNR after notch
notch_error = notch_filtered - clean;
snr_after_notch = 10*log10(sum(clean.^2)/sum(notch_error.^2));

%% ====== WIENER FILTER ======
% Calculate FFTs
N = fft(noise);
Y = fft(notch_filtered);

% Wiener filter
Pnn = abs(N).^2;
Pyy = abs(Y).^2;
H = max(0, 1 - Pnn./(Pyy + eps));

% Apply filter
filtered_freq = H .* Y;
wiener_filtered = real(ifft(filtered_freq));

% Final SNR
wiener_error = wiener_filtered - clean;
snr_after_wiener = 10*log10(sum(clean.^2)/sum(wiener_error.^2));

%% ====== PLOTTING ======
figure('Position', [100, 100, 1200, 800]);

% Original Signal
subplot(3,2,1);
plot(t, noisy);
title('Time: Noisy Signal');
xlabel('Time (s)');

subplot(3,2,2);
Y_orig = fft(noisy);
plot(f, abs(Y_orig(1:L/2+1))); % Fixed indexing
title('FFT: Noisy Signal');
xlabel('Frequency (Hz)');

% Notch Filtered
subplot(3,2,3);
plot(t, notch_filtered);
title('Time: After Notch Filtering');
xlabel('Time (s)');

subplot(3,2,4);
Y_notch = fft(notch_filtered);
plot(f, abs(Y_notch(1:L/2+1))); % Fixed indexing
title('FFT: After Notch Filtering');
xlabel('Frequency (Hz)');

% Wiener Filtered
subplot(3,2,5);
plot(t, wiener_filtered);
title('Time: After Wiener Filter');
xlabel('Time (s)');

subplot(3,2,6);
Y_wiener = fft(wiener_filtered);
plot(f, abs(Y_wiener(1:L/2+1))); % Fixed indexing
title('FFT: After Wiener Filter');
xlabel('Frequency (Hz)');

%% ====== RESULTS ======
fprintf('\nSNR Improvement:\n');
fprintf('1. Original: %.2f dB\n', snr_before);
fprintf('2. After Notch: %.2f dB (Δ%.1f dB)\n', snr_after_notch, snr_after_notch-snr_before);
fprintf('3. After Wiener: %.2f dB (Δ%.1f dB)\n', snr_after_wiener, snr_after_wiener-snr_after_notch);
