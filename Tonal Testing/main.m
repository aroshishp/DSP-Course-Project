clear; clc; close all;

noisy_speech = load('noisy_speech.txt'); 
external_noise = load('external_noise.txt'); 
clean_speech = load('clean_speech.txt');

fs = 44100; % Sampling frequency in Hz
t = 1/fs:1/fs:length(noisy_speech)/fs; % Time vector (2 seconds duration)

%% USER SETTINGS
added_freq = 2000;     % New programmable tone frequency
added_amp = 0.3;       % Amplitude of the added tone
Q_new = 60;            % Quality factor for new notch

% Add tone
[noisy_speech, external_noise] = add_tonal_noise(noisy_speech, external_noise, added_freq, added_amp, fs);


%% original SNR
signal_power = sum(clean_speech.^2);
noise_speech_power = sum((clean_speech - noisy_speech).^2);
Original_SNR = 10 * log10(signal_power / noise_speech_power);
fprintf('Original SNR: %.2f dB\n', Original_SNR);

% Parameters


%% RLS, Full-Suppression Mode
fs = 44100;
num_taps = 8; % Filter order
lambda = 0.99999; % Forgetting factor

fprintf('Two-Stage Full Suppression\n');

% First, apply only NEW frequency notch to external_noise
[filtered_noise_1, filtered_noise_2] = split_notch_filter(external_noise, fs, added_freq, Q_new);

% intermediate_signal = rls_filter(noisy_speech, filtered_noise_2, num_taps, lambda, 0);
% final_signal = rls_filter(intermediate_signal, filtered_noise_1, num_taps, lambda, 0);

% Stage 1: RLS with filtered_noise_1
% intermediate_signal = rls_filter(noisy_speech, filtered_noise_1, num_taps, 1, 0);

arr = 0.3 * sin(2 * pi * 2000 * t);
intermediate_signal = rls_filter(noisy_speech, arr', num_taps, 1, 0);
% Stage 2: RLS with filtered_noise_2
% final_signal = rls_filter(intermediate_signal, filtered_noise_2, num_taps, 1, 0);
final_signal = rls_filter(intermediate_signal, filtered_noise_1 , num_taps, 1, 0);

% Save and analyze
plotFFT(final_signal, fs);
plot_spectrogram(intermediate_signal, fs, 2)
plot_spectrogram(final_signal, fs, 2);
audiowrite('full_dual_rls.wav', final_signal, fs);
fprintf('Saved as full_dual_rls.wav\n');

% SNR
snr_full = 10 * log10(sum(clean_speech.^2) / sum((clean_speech - final_signal).^2));
fprintf('SNR after Full Dual-RLS: %.2f dB\n', snr_full);

% 
% % Save Audio
% audiowrite('rls_output.wav', enhanced_signal, fs);
% fprintf('Saved as rls_output.wav\n');


fprintf('Dual Notch Partial Suppression\n');
% Apply notch filters for both 1000 Hz and added_freq
temp = practical_notch_filter(external_noise, fs, 1000, 35);
ref_partial = practical_notch_filter(temp, fs, added_freq, Q_new);

% Run RLS once
partial_output = rls_filter(noisy_speech, ref_partial, num_taps, lambda, 0);

% Normalize and save
partial_output = partial_output / max(abs(partial_output));
audiowrite('partial_dual_notch.wav', partial_output, fs);
plotFFT(partial_output, fs);
plot_spectrogram(partial_output, fs, 2);
fprintf('Saved as partial_dual_notch.wav\n');
