clear all;
close all;

% Load the data files
fileNames = {'clean_speech.txt', 'external_noise.txt', 'noisy_speech.txt'};
clean_speech = load(fileNames{1});
noise_ref = load(fileNames{2});
noisy_speech = load(fileNames{3});

% Check if data has multiple columns, use only the first column
if size(clean_speech, 2) > 1
    clean_speech = clean_speech(:,1);
end
if size(noise_ref, 2) > 1
    noise_ref = noise_ref(:,1);
end
if size(noisy_speech, 2) > 1
    noisy_speech = noisy_speech(:,1);
end

% Ensure all signals have the same length
min_length = min([length(clean_speech), length(noise_ref), length(noisy_speech)]);
clean_speech = clean_speech(1:min_length);
noise_ref = noise_ref(1:min_length);
noisy_speech = noisy_speech(1:min_length);

% Parameters
fs = 44100;              % Sampling rate
filter_length = 512;      % Adaptive filter length
mu = 0.01;               % Initial step size for adaptive LMS
rho = 0.99;              % Forgetting factor for adaptive notch filter
num_samples = length(noisy_speech);

% Compute Initial SNR
snr_original = 10*log10(sum(clean_speech.^2) / sum((noisy_speech - clean_speech).^2));

% Adaptive Notch Filter Implementation
notch_freq = 1000; % Target frequency to remove
w0 = 2 * pi * notch_freq / fs;
B = [1, -2*cos(w0), 1]; % Notch filter coefficients
A = [1, -2*rho*cos(w0), rho^2];
notch_output = filter(B, A, noisy_speech);

% Compute SNR after Notch Filtering
snr_notch = 10*log10(sum(clean_speech.^2) / sum((notch_output - clean_speech).^2));

%  NLMS Filtering Implementation
w = zeros(filter_length, 1);  % Adaptive filter weights
e = zeros(num_samples, 1);    % Error signal
y = zeros(num_samples, 1);
epsilon = 1e-6;

for n = filter_length:num_samples
    x = noise_ref(n:-1:n-filter_length+1);
    y(n) = w' * x;
    e(n) = notch_output(n) - y(n);
    norm_factor = (x' * x) + epsilon;
    w = w + (2 * mu * e(n) * x) / norm_factor;
end

% Compute SNR after NLMS Filtering
snr_nlms = 10*log10(sum(clean_speech.^2) / sum((e - clean_speech).^2));

% Normalize signals
e = e / max(abs(e));

% Save as WAV files
audiowrite('cleaned_speech_nlms_notch.wav', e, fs);
audiowrite('notch_output.wav', notch_output, fs);

% Print SNR Report

fprintf('1. Original Noisy Speech SNR: %.2f dB\n', snr_original);
fprintf('2. After Adaptive Notch Filtering SNR: %.2f dB (Improvement: %.2f dB)\n', snr_notch, snr_notch - snr_original);
fprintf('3. After  NLMS Filtering SNR: %.2f dB (Improvement: %.2f dB)\n', snr_nlms, snr_nlms - snr_notch);
fprintf('4. Total Improvement: %.2f dB\n', snr_nlms - snr_original);

% Plot Signals and Spectrums
t = (0:num_samples-1) / fs;
f = linspace(0, fs/2, floor(num_samples/2)+1);

figure;
subplot(3,2,1); plot(t, noisy_speech); title('Noisy Speech'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,2,2); fft_noisy = fft(noisy_speech); plot(f, abs(fft_noisy(1:length(f)))./num_samples); title('Noisy Speech Spectrum'); xlabel('Frequency (Hz)');
subplot(3,2,3); plot(t, notch_output); title('After Adaptive Notch Filtering'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,2,4); fft_notch = fft(notch_output); plot(f, abs(fft_notch(1:length(f)))./num_samples); title('Spectrum After Notch'); xlabel('Frequency (Hz)');
subplot(3,2,5); plot(t, e); title('After  NLMS Filtering'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,2,6); fft_cleaned = fft(e); plot(f, abs(fft_cleaned(1:length(f)))./num_samples); title('Cleaned Speech Spectrum After NLMS'); xlabel('Frequency (Hz)');

% Play cleaned speech
fprintf('Playing cleaned speech after Adaptive Notch + NLMS...\n');
sound(e, fs);