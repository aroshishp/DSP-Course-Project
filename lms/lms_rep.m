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
filter_length = 64;      % Length of the adaptive filter (increased for better performance)
mu = 0.001;              % Step size (learning rate - reduced for stability)
num_samples = length(noisy_speech);

% Initialize
w = zeros(filter_length, 1);  % Filter weights
y = zeros(num_samples, 1);    % Filter output
e1 = zeros(num_samples, 1);    % Error signal (cleaned speech - 1st filtering)
e2 = zeros(num_samples, 1);    % Error signal (cleaned speech - 2nd filtering)


% LMS Algorithm - First Filtering
for n = filter_length:num_samples
    % Extract a window of the noise reference
    x = noise_ref(n:-1:n-filter_length+1);
    
    % Calculate filter output
    y(n) = w' * x;
    
    % Calculate error (this is our cleaned signal)
    e1(n) = noisy_speech(n) - y(n);
    
    % Update weights
    w = w + 2 * mu * e1(n) * x;
end

% Replace the beginning part that couldn't be processed due to filter length
e1(1:filter_length-1) = noisy_speech(1:filter_length-1);


% LMS Algorithm - Second Filtering (Applying LMS again on the first filtered signal)
w = zeros(filter_length, 1); % Re-initialize filter weights
y = zeros(num_samples, 1);    % Reset filter output
for n = filter_length:num_samples
    % Extract a window of the noise reference
    x = noise_ref(n:-1:n-filter_length+1);
    
    % Calculate filter output
    y(n) = w' * x;
    
    % Calculate error (this is our cleaned signal)
    e2(n) = e1(n) - y(n);  % Filtering the already filtered signal (e1)
    
    % Update weights
    w = w + 2 * mu * e2(n) * x;
end

% Replace the beginning part that couldn't be processed due to filter length
e2(1:filter_length-1) = e1(1:filter_length-1);


% Normalize signals to avoid clipping during playback and WAV file creation
clean_speech_norm = clean_speech / max(abs(clean_speech));
noisy_speech_norm = noisy_speech / max(abs(noisy_speech));
cleaned_speech_norm1 = e1 / max(abs(e1)); % Normalize the first filtered signal
cleaned_speech_norm2 = e2 / max(abs(e2)); % Normalize the second filtered signal

% Save as WAV files
audiowrite('clean_speech.wav', clean_speech_norm, fs);
audiowrite('noisy_speech.wav', noisy_speech_norm, fs);
audiowrite('cleaned_speech_1.wav', cleaned_speech_norm1, fs); % Save the first filtered output
audiowrite('cleaned_speech_2.wav', cleaned_speech_norm2, fs); % Save the second filtered output


fprintf('\nAudio files saved as WAV files:\n');
fprintf('1. clean_speech.wav - Original clean speech\n');
fprintf('2. noisy_speech.wav - Noisy speech\n');
fprintf('3. cleaned_speech_1.wav - Cleaned speech after *one* filtering\n');
fprintf('4. cleaned_speech_2.wav - Cleaned speech after *two* filtering passes\n');

% Time vectors for plotting
t = (0:num_samples-1) / fs;

% Compute Fourier Transforms
N = length(clean_speech);
f = (0:N-1) * (fs / N); % Frequency axis
fft_clean = abs(fft(clean_speech));
fft_noisy = abs(fft(noisy_speech));
fft_cleaned1 = abs(fft(e1)); % FFT of 1st filtered signal
fft_cleaned2 = abs(fft(e2)); % FFT of 2nd filtered signal


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
plot(t, cleaned_speech_norm1);
title('Cleaned Speech Waveform (1 filter pass)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,2,6);
plot(f(1:N/2), fft_cleaned1(1:N/2));
title('Cleaned Speech Spectrum (1 filter pass)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

figure;  % create new figure for the 2nd filter
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
plot(t, cleaned_speech_norm2);
title('Cleaned Speech Waveform (2 filter passes)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,2,6);
plot(f(1:N/2), fft_cleaned2(1:N/2));
title('Cleaned Speech Spectrum (2 filter passes)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


% Calculate SNR improvement (for first pass)
noisy_snr = 10*log10(sum(clean_speech.^2)/sum((noisy_speech-clean_speech).^2));
cleaned_snr1 = 10*log10(sum(clean_speech.^2)/sum((e1-clean_speech).^2));
improvement1 = cleaned_snr1 - noisy_snr;

% Calculate SNR improvement (for second pass)
cleaned_snr2 = 10*log10(sum(clean_speech.^2)/sum((e2-clean_speech).^2));
improvement2 = cleaned_snr2 - noisy_snr;


fprintf('SNR before noise cancellation: %.2f dB\n', noisy_snr);
fprintf('SNR after *one* noise cancellation pass: %.2f dB\n', cleaned_snr1);
fprintf('SNR improvement (1 pass): %.2f dB\n', improvement1);
fprintf('SNR after *two* noise cancellation passes: %.2f dB\n', cleaned_snr2);
fprintf('SNR improvement (2 passes): %.2f dB\n', improvement2);


% Play audio in sequence with proper pauses
fprintf('\nPlaying original clean speech...\n');
%sound(clean_speech_norm, fs);
%pause(length(clean_speech_norm)/fs + 1);

fprintf('Playing noisy speech...\n');
%sound(noisy_speech_norm, fs);
%pause(length(noisy_speech_norm)/fs + 1);

fprintf('Playing cleaned speech (1 pass)...\n');
sound(cleaned_speech_norm1, fs);
pause(length(cleaned_speech_norm1)/fs + 1);


fprintf('Playing cleaned speech (2 passes)...\n');
sound(cleaned_speech_norm2, fs);
