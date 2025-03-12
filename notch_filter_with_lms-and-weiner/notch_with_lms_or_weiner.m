%% Load Signals
clean = load('clean_speech.txt');
noise = load('external_noise.txt');
noisy = load('noisy_speech.txt');

fs = 44100; 
t = (0:length(clean)-1)/fs;

%% Signal Alignment and Normalization
min_len = min([length(clean), length(noise), length(noisy)]);
clean = clean(1:min_len);
noise = noise(1:min_len);
noisy = noisy(1:min_len);

clean = clean/max(abs(clean));
noise = noise/max(abs(noise));
noisy = noisy/max(abs(noisy));

%% Initialize Figure
figure('Name','Enhanced Processing Pipeline', 'Position', [100 100 1400 1200])

%% Stage 0: Clean Signal Reference
subplot(5,2,1)
plot_fft(clean, fs, gca, 'Clean Reference Spectrum', 'b', '-')
subplot(5,2,2)
plot(t, clean)
title('Time Domain - Clean Speech')
xlabel('Time (s)'), axis tight

%% Stage 1: Noisy Signal Analysis
subplot(5,2,3)
plot_fft(noisy, fs, gca, 'Original Noisy Spectrum', 'r', '-')
subplot(5,2,4)
plot(t, noisy)
title('Time Domain - Noisy Speech')
xlabel('Time (s)'), axis tight

%% Improved Notch Filtering
Y_noisy = fft(noisy);
L = length(noisy);
f = fs*(0:(floor(L/2)))/L;
spectrum = abs(Y_noisy(1:floor(L/2)+1));

threshold = 0.2 * max(spectrum);
[pks,locs] = findpeaks(spectrum, 'MinPeakHeight', threshold, 'MinPeakDistance', 50);

tonal_freqs = f(locs);
valid_freqs = tonal_freqs(tonal_freqs < fs/2*0.99);

filtered = noisy;
for k = 1:length(valid_freqs)
    wo = valid_freqs(k)/(fs/2);
    bw = 75/(fs/2);
    [b,a] = iirnotch(wo, bw);
    filtered = filter(b,a,filtered);
end

subplot(5,2,5)
plot_fft(filtered, fs, gca, 'After Notch Filtering', 'm', '-')
hold on
plot_fft(clean, fs, gca, 'Clean Reference', 'b', ':')  % Fixed line
hold off
subplot(5,2,6)
plot(t, filtered)
title('Time Domain - Notched Speech')
xlabel('Time (s)'), axis tight

%% Stage 2: Wiener Filtering
H_wiener = (abs(fft(clean)).^2) ./ (abs(fft(clean)).^2 + abs(fft(noise)).^2);
wiener_filtered = real(ifft(H_wiener .* fft(filtered)));

subplot(5,2,7)
plot_fft(wiener_filtered, fs, gca, 'After Wiener Filtering', 'g', '-')
hold on
plot_fft(clean, fs, gca, 'Clean Reference', 'b', ':')  % Fixed line
hold off
subplot(5,2,8)
plot(t, wiener_filtered)
title('Time Domain - Wiener Filtered')
xlabel('Time (s)'), axis tight

%% Stage 3: LMS Adaptive Filter
order = 128;
mu = 0.001;
lms = dsp.LMSFilter(order, 'StepSize', mu, 'Method', 'Normalized LMS');
[~, lms_filtered] = lms(noise(1:length(wiener_filtered)), wiener_filtered);

subplot(5,2,9)
plot_fft(clean, fs, gca, 'Clean Reference', 'b', ':')  % Fixed line
hold on
plot_fft(lms_filtered, fs, gca, 'After LMS Filtering', 'r', '-')
hold off
subplot(5,2,10)
plot(t, lms_filtered)
title('Time Domain - LMS Filtered')
xlabel('Time (s)'), axis tight

%% SNR Calculations
snr_results = zeros(5,1);
snr_results(1) = calculate_snr(clean, noisy);
snr_results(2) = calculate_snr(clean, filtered);
snr_results(3) = calculate_snr(clean, wiener_filtered);
snr_results(4) = calculate_snr(clean, lms_filtered);

fprintf('\nSNR Improvement Report:\n');
fprintf('----------------------------------\n');
fprintf('| Processing Stage   | SNR (dB) |\n');
fprintf('|--------------------|----------|\n');
fprintf('| Clean Reference    |   ∞      |\n');
fprintf('| Original Noisy     | %7.2f  |\n', snr_results(1));
fprintf('| After Notch        | %7.2f  |\n', snr_results(2));
fprintf('| After Wiener       | %7.2f  |\n', snr_results(3));
fprintf('| After LMS          | %7.2f  |\n', snr_results(4));
fprintf('----------------------------------\n');

%% Audio Playback
fprintf('\nPlaying processing stages...\n');
soundsc(clean, fs); pause(length(clean)/fs + 1)
soundsc(noisy, fs); pause(length(noisy)/fs + 1)
soundsc(filtered, fs); pause(length(filtered)/fs + 1)
soundsc(wiener_filtered, fs); pause(length(wiener_filtered)/fs + 1)
soundsc(lms_filtered, fs);

%% FFT Analysis Function
function plot_fft(signal, fs, ax, title_str, color, linestyle)
    if nargin < 6
        linestyle = '-';
    end
    L = length(signal);
    f = fs*(0:(floor(L/2)))/L;
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    plot(ax, f, P1, 'Color', color, 'LineStyle', linestyle)
    title(ax, title_str)
    xlabel(ax, 'Frequency (Hz)')
    ylabel(ax, 'Magnitude')
    xlim([0 fs/2])
    grid(ax, 'on')
end

%% SNR Calculation Function
function snr = calculate_snr(clean_sig, processed_sig)
    noise_component = clean_sig - processed_sig;
    signal_power = sum(clean_sig.^2);
    noise_power = sum(noise_component.^2);
    snr = 10*log10(signal_power/noise_power);
end
