noisy_speech = load('noisy_speech.txt');
external_noise = load('external_noise.txt');
clean_speech = load('clean_speech.txt');

% Parameters
fs = 44100; 
tones = [1000]; 
num_taps = 8; 

signal_power = sum(clean_speech.^2);
noise_speech_power = sum((clean_speech - noisy_speech).^2);
Original_SNR = 10 * log10(signal_power / noise_speech_power);

lambda_values = [0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 0.9998, 0.9999, 0.99995, 0.99998, 0.99999, 0.999995];
num_lambdas = length(lambda_values);

snr_results = zeros(num_lambdas, 1);

for i = 1:num_lambdas
    lambda = lambda_values(i);
    
    enhanced_signal = rls_filter(noisy_speech, external_noise, fs, num_taps, lambda, Q, 0, tones);

    rls_noise_power = sum((clean_speech - enhanced_signal).^2);
    snr_results(i) = 10 * log10(signal_power / rls_noise_power);
end

figure;
plot(lambda_values, snr_results, '-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Lambda');
ylabel('SNR (dB)');
xlim([0.9, 1]);
title('SNR vs Forgetting Factor');


