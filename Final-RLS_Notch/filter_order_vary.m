noisy_speech = load('noisy_speech.txt');
external_noise = load('external_noise.txt');
clean_speech = load('clean_speech.txt');

% Parameters
fs = 44100; 
tones = [1000]; 
lambda = 0.99999; 

filter_orders = unique(round(logspace(log10(2), log10(128), 15)));
num_orders = length(filter_orders);

snr_results = zeros(num_orders, 1);
compute_time = zeros(num_orders, 1);

for i = 1:num_orders
    num_taps = filter_orders(i);
   
    tic;
    enhanced_signal = rls_filter(noisy_speech, external_noise, fs, num_taps, lambda, Q, 0, tones);
    compute_time(i) = toc;

    rls_noise_power = sum((clean_speech - enhanced_signal).^2);
    snr_results(i) = 10 * log10(signal_power / rls_noise_power);
end

figure;
plot(filter_orders, snr_results, '-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Filter Order');
ylabel('SNR (dB)');
title('RLS Performance vs Filter Order');
hold on;
plot([min(filter_orders), max(filter_orders)], [Original_SNR, Original_SNR], '--r', 'LineWidth', 1.5);
legend('RLS SNR', 'Original SNR');

% Plot computation time
figure;
plot(filter_orders, compute_time, '-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Filter Order');
ylabel('Computation Time (s)');
title('RLS Computation Time vs Filter Order');
