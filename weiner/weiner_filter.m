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
C = fft(clean);
N = fft(noise);
Y = fft(noisy);

% Calculate Wiener filter
H = (abs(C).^2) ./ (abs(C).^2 + abs(N).^2);
filtered_freq = H .* Y;
filtered = real(ifft(filtered_freq));

% Play signals (adjust fs according to your sampling rate)
fs = 44100; % Typical sampling rate for speech
soundsc(clean, fs); 
pause(length(clean)/fs + 1);
soundsc(filtered, fs);

% Calculate SNR metrics
original_noise = noisy - clean;
snr_before = 10*log10(sum(clean.^2)/sum(original_noise.^2));
residual_noise = clean - filtered;
snr_after = 10*log10(sum(clean.^2)/sum(residual_noise.^2));

fprintf('SNR before filtering: %.2f dB\n', snr_before);
fprintf('SNR after filtering: %.2f dB\n', snr_after);
