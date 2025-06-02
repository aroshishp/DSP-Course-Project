noisy_speech = load('noisy_speech.txt');
external_noise = load('external_noise.txt');
clean_speech = load('clean_speech.txt');

% Parameters
fs = 44100; 
tones = [1000]; 
num_taps = 8; 
lambda = 0.99999; 

Q_values = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500];
num_Q = length(Q_values);
   
tpr_results = zeros(num_Q, 1);    
spr_results = zeros(num_Q, 1);    

for i = 1:num_Q
    Q = Q_values(i);

    partial_notched = rls_filter(noisy_speech, external_noise, fs, num_taps, lambda, Q, 1, tones);

    tpr_results(i) = tone_preservation_ratio(noisy_speech, partial_notched, 1000, fs);
    spr_results(i) = spectrum_preservation_ratio(clean_speech, partial_notched, 1000, fs);
end

figure;
yyaxis left;
plot(Q_values, tpr_results, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Tone Preservation Index (TPR)');
yyaxis right;
plot(Q_values, spr_results, '-s', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Spectrum Preservation Ratio (SPR)');
grid on;
xlabel('Quality Factor (Q)');
title('Preservation Metrics vs Quality Factor - Partial Suppression');
legend('TPR', 'SPR');