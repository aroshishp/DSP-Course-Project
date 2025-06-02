noisy_speech = load('noisy_speech.txt'); 
external_noise = load('external_noise.txt'); 
clean_speech = load('clean_speech.txt');

signal_power = sum(clean_speech.^2);
noise_speech_power = sum((clean_speech - noisy_speech).^2);
Original_SNR = 10 * log10(signal_power / noise_speech_power);
fprintf('Original SNR: %.2f dB\n', Original_SNR);

fs = 44100;

num_taps = 2;
lambda = 0.99999; 
Q = 35;
mu = 0.001;              
num_samples = length(noisy_speech);
epsilon = 1e-10;

w_nlms = zeros(num_taps, 1);      
enhanced_signal = zeros(size(noisy_speech)); 
filtered_buffer = zeros(size(external_noise)); 

for n = 1:length(noisy_speech)
    filtered_buffer(n) = external_noise(n); 

    if n >= num_taps
        x_vec = filtered_buffer(n - num_taps + 1 : n);
        x_vec = flip(x_vec);
        d = noisy_speech(n);
        y_hat = w_nlms' * x_vec;
        e = d - y_hat;
        norm_factor = (x_vec' * x_vec) + epsilon;
        w_nlms = w_nlms + (mu / norm_factor) * x_vec * e;
        enhanced_signal(n) = e;
    end
end

rls_noise_power = sum((clean_speech - enhanced_signal).^2);
SNR_RLS = 10 * log10(signal_power / rls_noise_power);
fprintf('SNR after NLMS: %.2f dB\n', SNR_RLS);