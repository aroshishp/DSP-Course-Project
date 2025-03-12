clear; clc; close all;

%% Load Input Signals
noisy_speech = load('noisy_speech.txt'); % s(n) + v(n)
external_noise = load('external_noise.txt'); % w(n)
clean_speech = load('clean_speech.txt');

fs = 44100;  % Sampling rate 
frame_size = 2048; % Buffer size for real-time processing
num_taps = 64; % Filter order for adaptive RLS
lambda = 0.99; % Forgetting factor for RLS

%% Initialize RLS Adaptive Filter
delta = 0.01; % Small positive constant to initialize P matrix
w_rls = zeros(num_taps, 1); % Filter coefficients
P = (1/delta) * eye(num_taps); % Inverse correlation matrix

%% Notch Filter Design (1000 Hz)
wo = 1000 / (fs/2); % Normalized frequency
bw = wo / 35; % Bandwidth of notch
[b_notch, a_notch] = iirnotch(wo, bw);

%% Output Storage
enhanced_signal = zeros(size(noisy_speech));

%% Real-Time Adaptive Filtering (Frame-by-Frame Processing)
for i = 1:frame_size:length(noisy_speech)-frame_size
    % Get current frame
    x = external_noise(i:i+frame_size-1); % Reference noise w(n)
    d = noisy_speech(i:i+frame_size-1);  % Noisy speech s(n) + v(n)

    % Construct Input Matrix for Adaptive Filtering
    X = zeros(frame_size, num_taps);
    for j = num_taps:frame_size
        X(j, :) = flip(x(j-num_taps+1:j)); % Delayed reference noise samples
    end
    
    % RLS Adaptive Filtering
    for j = num_taps:frame_size
        x_j = X(j, :)'; % Current reference signal segment
        K = (P * x_j) / (lambda + x_j' * P * x_j); % Gain vector
        v_hat = w_rls' * x_j; % Estimated noise v̂(n)
        e = d(j) - v_hat; % Error signal (speech estimate)
        w_rls = w_rls + K * e; % Update filter weights
        P = (P - K * x_j' * P) / lambda; % Update inverse correlation matrix
        enhanced_signal(i + j - 1) = e; % Store enhanced speech
    end
end

%% Apply Notch Filter (Removes Tonal Noise at 1000 Hz)
enhanced_signal = filter(b_notch, a_notch, enhanced_signal);

%% Apply Spectral Subtraction
frame_overlap = frame_size / 2;
hann_win = hann(frame_size);

enhanced_spectral = zeros(size(enhanced_signal));
for i = 1:frame_overlap:length(enhanced_signal)-frame_size
    frame = enhanced_signal(i:i+frame_size-1) .* hann_win;
    Y = fft(frame);

    % Estimate Noise Spectrum (Assume first few frames are noise only)
    if i < 10 * frame_overlap
        noise_mag = abs(Y); % Noise magnitude estimate
    end

    % Spectral Subtraction
    Y_clean = Y - 1.2 * noise_mag; % Subtract noise estimate
    Y_clean(Y_clean < 0) = 0; % Prevent negative values

    % Reconstruct Signal
    frame_clean = real(ifft(Y_clean)) .* hann_win;
    enhanced_spectral(i:i+frame_size-1) = enhanced_spectral(i:i+frame_size-1) + frame_clean;
end


%% Calculate SNR
signal_power = sum(clean_speech.^2);
noise_power = sum((clean_speech - enhanced_spectral).^2);
SNR = 10 * log10(signal_power / noise_power);
fprintf('Final SNR: %.2f dB\n', SNR);


%% Save and Plot Results
audiowrite('enhanced_output.wav', enhanced_spectral, fs);
figure();
plot(enhanced_spectral);
title('Enhanced Speech Signal');
xlabel('Sample Index'); ylabel('Amplitude');
legend('Processed Signal');

figure();
plot(clean_speech);
