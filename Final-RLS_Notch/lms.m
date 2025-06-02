noisy_speech = load('noisy_speech.txt'); 
external_noise = load('external_noise.txt'); 
clean_speech = load('clean_speech.txt');

% Parameters
fs = 44100;              
filter_length = 8;      
mu = 0.001;              
num_samples = length(noisy_speech);

w = zeros(filter_length, 1);  
y = zeros(num_samples, 1);    
e = zeros(num_samples, 1);    
alpha = 0.01;

for n = filter_length:num_samples
    x = external_noise(n:-1:n-filter_length+1);
    
    y(n) = w' * x;
    e(n) = noisy_speech(n) - y(n);
    
    w = w + 2 * mu * e(n) * x;
end

clean_speech_norm = clean_speech / max(abs(clean_speech));
noisy_speech_norm = noisy_speech / max(abs(noisy_speech));
Filtered_speech_norm = e / max(abs(e));

noisy_snr = 10*log10(sum(clean_speech.^2)/sum((noisy_speech-clean_speech).^2));
cleaned_snr = 10*log10(sum(clean_speech.^2)/sum((e-clean_speech).^2));
improvement = cleaned_snr - noisy_snr;

fprintf('SNR before noise cancellation: %.2f dB\n', noisy_snr);
fprintf('SNR after noise cancellation: %.2f dB\n', cleaned_snr);
fprintf('SNR improvement: %.2f dB\n', improvement);


