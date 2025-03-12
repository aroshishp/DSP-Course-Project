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
fs = 44100;              % Sampling rate (matching your code)
filter_length = 512;      % Length of the adaptive filter (increased for better performance)
mu = 0.0005;              % Step size (learning rate - reduced for stability)
num_samples = length(noisy_speech);

% Initialize
w = zeros(filter_length, 1);  % Filter weights
y = zeros(num_samples, 1);    % Filter output
e = zeros(num_samples, 1);    % Error signal (cleaned speech)

alpha = 0.001;
% LMS Algorithm
for n = filter_length:num_samples
    % Extract a window of the noise reference
    x = noise_ref(n:-1:n-filter_length+1);
    
    % Calculate filter output
    y(n) = w' * x;
    
    % Calculate error (this is our cleaned signal)
    e(n) = noisy_speech(n) - y(n);
    
    % Update weights
    % w = w + 2 * mu * e(n) * x;
    w = (1 - alpha*mu) * w + 2 * mu * e(n) * x; % alpha is small constant.
end

% Replace the beginning part that couldn't be processed due to filter length
e(1:filter_length-1) = noisy_speech(1:filter_length-1);

% Normalize signals to avoid clipping during playback and WAV file creation
clean_speech_norm = clean_speech / max(abs(clean_speech));
noisy_speech_norm = noisy_speech / max(abs(noisy_speech));
cleaned_speech_norm = e / max(abs(e));

% Save as WAV files
audiowrite('clean_speech.wav', clean_speech_norm, fs);
audiowrite('noisy_speech.wav', noisy_speech_norm, fs);
audiowrite('cleaned_speech.wav', cleaned_speech_norm, fs);

fprintf('\nAudio files saved as WAV files:\n');
fprintf('1. clean_speech.wav - Original clean speech\n');
fprintf('2. noisy_speech.wav - Noisy speech\n');
fprintf('3. cleaned_speech.wav - Cleaned speech after noise cancellation\n');

% Time vectors for plotting
t = (0:num_samples-1) / fs;

% Compute Fourier Transforms
N = length(clean_speech);
f = (0:N-1) * (fs / N); % Frequency axis
fft_clean = abs(fft(clean_speech));
fft_noisy = abs(fft(noisy_speech));
fft_cleaned = abs(fft(e));

% Plot waveforms and frequency spectra
figure;
subplot(3,2,1);
plot(t, clean_speech_norm);
title('Clean Speech Waveform');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,2,2);
plot(f(1:N/2), fft_clean(1:N/2));
title('Clean Speech Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,2,3);
plot(t, noisy_speech_norm);
title('Noisy Speech Waveform');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,2,4);
plot(f(1:N/2), fft_noisy(1:N/2));
title('Noisy Speech Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,2,5);
plot(t, cleaned_speech_norm);
title('Cleaned Speech Waveform');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,2,6);
plot(f(1:N/2), fft_cleaned(1:N/2));
title('Cleaned Speech Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Calculate SNR improvement
noisy_snr = 10*log10(sum(clean_speech.^2)/sum((noisy_speech-clean_speech).^2));
cleaned_snr = 10*log10(sum(clean_speech.^2)/sum((e-clean_speech).^2));
improvement = cleaned_snr - noisy_snr;

fprintf('SNR before noise cancellation: %.2f dB\n', noisy_snr);
fprintf('SNR after noise cancellation: %.2f dB\n', cleaned_snr);
fprintf('SNR improvement: %.2f dB\n', improvement);

% Play audio in sequence with proper pauses
fprintf('\nPlaying original clean speech...\n');
% sound(clean_speech_norm, fs);
% pause(length(clean_speech_norm)/fs + 1);

fprintf('Playing noisy speech...\n');
sound(noisy_speech_norm, fs);
pause(length(noisy_speech_norm)/fs + 1);

fprintf('Playing cleaned speech...\n');
sound(cleaned_speech_norm, fs);
