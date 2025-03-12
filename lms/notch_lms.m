clear all;
close all;

% Load data files
files = {'clean_speech.txt', 'external_noise.txt', 'noisy_speech.txt'};
clean = load(files{1});
noise_ref = load(files{2});
noisy = load(files{3});

% Ensure single-channel audio
clean = clean(:,1);
noise_ref = noise_ref(:,1);
noisy = noisy(:,1);

% Match lengths
L = min([length(clean), length(noise_ref), length(noisy)]);
clean = clean(1:L);
noise_ref = noise_ref(1:L);
noisy = noisy(1:L);

%% ====== ADAPTIVE NOTCH FILTERING ======
fs = 44100;  % Sample rate

% Frequency analysis for tone detection
Y_noisy = fft(noisy);
f = fs*(0:(floor(L/2)))/L;
spectrum = abs(Y_noisy(1:floor(L/2)+1));

% Peak detection parameters
threshold = 0.2 * max(spectrum);
min_peak_distance = 50;  % Minimum 50 Hz between peaks

% Find tonal components
[pks,locs] = findpeaks(spectrum,...
    'MinPeakHeight', threshold,...
    'MinPeakDistance', min_peak_distance);

% Filter parameters
bw = 25;  % Bandwidth in Hz
filtered = noisy;

% Apply multiple notch filters
for k = 1:length(locs)
    freq = f(locs(k));
    if freq < 0.99*(fs/2)  % Avoid Nyquist frequency
        wo = freq/(fs/2);
        bw_norm = bw/(fs/2);
        [b,a] = iirnotch(wo, bw_norm);
        filtered = filter(b,a,filtered);
    end
end

noisy_speech = filtered;

%% ====== LMS NOISE CANCELLATION ======
% Adaptive filter parameters
filter_order = 512;
mu = 0.0005;
alpha = 0.001;  % Leakage factor

% Initialize variables
w = zeros(filter_order, 1);
e = zeros(L,1);

% LMS algorithm
for n = filter_order:L
    x = noise_ref(n:-1:n-filter_order+1);
    y = w' * x;
    e(n) = noisy_speech(n) - y;
    w = (1 - alpha*mu)*w + 2*mu*e(n)*x;
end

% Replace unprocessed samples
e(1:filter_order-1) = noisy_speech(1:filter_order-1);

%% ====== ANALYSIS & OUTPUT ======
% Normalization
clean_norm = clean/max(abs(clean));
noisy_norm = noisy/max(abs(noisy));
notched_norm = noisy_speech/max(abs(noisy_speech));
cleaned_norm = e/max(abs(e));

% Save audio
audiowrite('0_clean.wav', clean_norm, fs);
audiowrite('1_noisy.wav', noisy_norm, fs);
audiowrite('2_notched.wav', notched_norm, fs);
audiowrite('3_cleaned.wav', cleaned_norm, fs);

% SNR calculations
original_snr = 10*log10(sum(clean.^2)/sum((noisy - clean).^2));
notch_snr = 10*log10(sum(clean.^2)/sum((noisy_speech - clean).^2));
final_snr = 10*log10(sum(clean.^2)/sum((e - clean).^2));

fprintf('SNR Improvement:\n');
fprintf('Original: %.2f dB\n', original_snr);
fprintf('After Notch: %.2f dB\n', notch_snr);
fprintf('Final Cleaned: %.2f dB\n', final_snr);
fprintf('Total Improvement: %.2f dB\n', final_snr - original_snr);

%% ====== VISUALIZATION ======
t = (0:L-1)/fs;

% Time domain plots
figure;
subplot(411); plot(t, clean); title('Clean Speech');
subplot(412); plot(t, noisy); title('Original Noisy');
subplot(413); plot(t, noisy_speech); title('Notch Filtered');
subplot(414); plot(t, e); title('Final Output');

% Frequency domain plots (CORRECTED VERSION)
figure;
hold on;
NFFT = 2^nextpow2(L);
f_axis = fs/2*linspace(0,1,NFFT/2+1);  % Renamed to avoid conflict with frequency variable

% Compute FFTs first
Y_clean = fft(clean, NFFT);
Y_noisy = fft(noisy, NFFT);
Y_notched = fft(noisy_speech, NFFT);
Y_cleaned = fft(e, NFFT);

% Plot using computed FFT results
plot(f_axis, abs(Y_clean(1:NFFT/2+1)), 'g');
% plot(f_axis, abs(Y_noisy(1:NFFT/2+1)), 'r');
% plot(f_axis, abs(Y_notched(1:NFFT/2+1)), 'm');
plot(f_axis, abs(Y_cleaned(1:NFFT/2+1)), 'b');

legend('Clean','Noisy','Notched','Cleaned');
xlabel('Frequency (Hz)');
ylabel('dB');
title('Frequency Spectrum Comparison');
grid on;
xlim([0 5000]);


% Play audio (uncomment if desired)
% sound(original_noisy_norm, fs);
% pause(length(original_noisy_norm)/fs + 1);
% sound(notch_filtered_norm, fs);
% pause(length(notch_filtered_norm)/fs + 1);
% sound(cleaned_speech_norm, fs);
