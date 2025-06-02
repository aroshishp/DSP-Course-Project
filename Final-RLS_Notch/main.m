clear; clc; close all;

noisy_speech = load('noisy_speech.txt'); 
external_noise = load('external_noise.txt'); 
clean_speech = load('clean_speech.txt');

%% original SNR
signal_power = sum(clean_speech.^2);
noise_speech_power = sum((clean_speech - noisy_speech).^2);
Original_SNR = 10 * log10(signal_power / noise_speech_power);
fprintf('Original SNR: %.2f dB\n\n', Original_SNR);

% Parameters
fs = 44100; %Sampling frequency

%% TONAL FREQUENCIES USER SUPPLIED
tones = [1000];

%% Full-Suppression Mode
num_taps = 8;       %Filter order
lambda = 0.99999;   %Forgetting factor
R = 0.999;          %Notch Filter Quality Factor for Partial Suppr.

fprintf('Full Suppression Mode\n');
enhanced_signal = rls_filter(noisy_speech, external_noise, fs, num_taps, lambda, R, 0, tones);

rls_noise_power = sum((clean_speech - enhanced_signal).^2);
SNR_RLS = 10 * log10(signal_power / rls_noise_power);

fprintf('SNR after Full Suppression: %.2f dB\n\n', SNR_RLS);

plot_spectrogram(clean_speech, fs, 0);
plot_spectrogram(enhanced_signal, fs, 1);

% % Save Audio
% enhanced_signal_norm = enhanced_signal / (max(abs(enhanced_signal)));
% audiowrite('rls_output.wav', enhanced_signal_norm, fs);
% fprintf('Saved as rls_output.wav\n');

%% Partial Suppression Mode
fprintf('Partial Suppression Mode\n');
partial_notched = rls_filter(noisy_speech, external_noise, fs, num_taps, lambda, R, 1, tones);

tpr = tone_preservation_ratio(noisy_speech, partial_notched, 1000, fs);
spr = spectrum_preservation_ratio(clean_speech, partial_notched, 1000, fs);

fprintf("Tone Preservation Ratio = %f\n", tpr);
fprintf("Non-Tone Spectrum Preservation Ratio = %f\n", spr);

plot_spectrogram(partial_notched, fs, 2);
plotFFT(partial_notched, fs, 1);

% % Save Audio
% partial_practical_notched = partial_practical_notched / (max(abs(partial_practical_notched)));
% audiowrite('practical_partial.wav', partial_practical_notched, fs);
% fprintf('Saved as practical_partial.wav\n');