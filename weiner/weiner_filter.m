%%%%% Weiner without spectral subtraction
% Load signals from files
clean = load('clean_speech.txt');
noise = load('external_noise.txt');
noisy = load('noisy_speech.txt');

% Ensure column vectors
clean = clean(:);
noise = noise(:);
noisy = noisy(:);

% Match signal lengths
min_len = min([length(clean), length(noise), length(noisy)]);
clean = clean(1:min_len);
noise = noise(1:min_len);
noisy = noisy(1:min_len);

% Compute FFTs
N = fft(noise);
Y = fft(noisy);

% Calculate Wiener filter using noise estimate and noisy signal
H = max(0, 1 - (abs(N).^2) ./ (abs(Y).^2 + eps));  % eps prevents division by zero
filtered_freq = H .* Y;
filtered = real(ifft(filtered_freq));

% Play signals (adjust fs according to your sampling rate)
fs = 44100; % Typical sampling rate for speech
% soundsc(noisy, fs); 
% pause(length(noisy)/fs + 1);
% soundsc(filtered, fs);

% Calculate SNR metrics (clean signal only used here)
original_noise = noisy - clean;
snr_before = 10*log10(sum(clean.^2)/sum(original_noise.^2));
residual_noise = clean - filtered;
snr_after = 10*log10(sum(clean.^2)/sum(residual_noise.^2));

fprintf('SNR before filtering: %.2f dB\n', snr_before);
fprintf('SNR after filtering without spectral subtraction: %.2f dB\n', snr_after);






%%%%% Weiner with spectral subtraction
% Ensure column vectors & match lengths
clean = clean(:); noise = noise(:); noisy = noisy(:);
min_len = min([length(clean), length(noise), length(noisy)]);
clean = clean(1:min_len); noise = noise(1:min_len); noisy = noisy(1:min_len);

% Compute PSD estimates
Y = fft(noisy);
N = fft(noise);

% Improved Wiener implementation using spectral subtraction
PSD_noise = abs(N).^2;
PSD_noisy = abs(Y).^2;
PSD_clean_estimate = max(PSD_noisy - PSD_noise, 0);  % Spectral subtraction
H = PSD_clean_estimate ./ (PSD_clean_estimate + PSD_noise + eps);

% Apply filtering
filtered_freq = H .* Y;
filtered = real(ifft(filtered_freq));

% Playback and SNR calculation (clean only used here)
fs = 44100;
% soundsc(noisy, fs); pause(length(noisy)/fs + 1);
% soundsc(filtered, fs);

% SNR metrics
original_noise = noisy - clean;
snr_before = 10*log10(sum(clean.^2)/sum(original_noise.^2));
residual_noise = clean - filtered;
snr_after = 10*log10(sum(clean.^2)/sum(residual_noise.^2));


fprintf('SNR after filtering with spectral subtraction: %.2f dB\n', snr_after);

