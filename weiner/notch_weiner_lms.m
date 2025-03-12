%% Load Signals
clean = load('clean_speech.txt');
noise = load('external_noise.txt');
noisy = load('noisy_speech.txt');

fs = 44100; 
t = (0:length(noisy)-1)/fs;

%% Signal Alignment
min_len = min([length(clean), length(noise), length(noisy)]);
clean = clean(1:min_len);
noise = noise(1:min_len);
noisy = noisy(1:min_len);

%% Initialize Figure
figure('Name','Enhanced Processing Pipeline', 'Position', [100 100 1400 900])

%% Stage 1: Noisy Signal Analysis
subplot(4,2,1)
plot_fft(noisy, fs, gca, 'Noisy Spectrum', 'r', '-')
subplot(4,2,2)
plot(t, noisy)
title('Time Domain - Noisy Speech')
xlabel('Time (s)'), axis tight

%% Improved Notch Filtering
[filtered, tonal_freqs] = adaptive_notch_filter(noisy, fs);

subplot(4,2,3)
plot_fft(filtered, fs, gca, 'After Notch', 'm', '-')
subplot(4,2,4)
plot(t, filtered)
title('Time Domain - Notched Speech')
xlabel('Time (s)'), axis tight

%% Stage 2: Wiener Filtering
H = wiener_filter(noisy, noise);
wiener_filtered = real(ifft(H .* fft(filtered)));

subplot(4,2,5)
plot_fft(wiener_filtered, fs, gca, 'After Wiener', 'g', '-')
subplot(4,2,6)
plot(t, wiener_filtered)
title('Time Domain - Wiener Filtered')
xlabel('Time (s)'), axis tight

%% Stage 3: LMS Adaptive Filter
lms_filtered = adaptive_lms(wiener_filtered, noise);

subplot(4,2,7)
plot_fft(lms_filtered, fs, gca, 'After LMS', 'b', '-')
subplot(4,2,8)
plot(t, lms_filtered)
title('Time Domain - LMS Filtered')
xlabel('Time (s)'), axis tight

%% SNR Calculations Using Clean Reference
snr_original = calculate_snr(clean, noisy);
snr_notch = calculate_snr(clean, filtered);
snr_wiener = calculate_snr(clean, wiener_filtered);
snr_lms = calculate_snr(clean, lms_filtered);

fprintf('\nSNR Improvement Report:\n');
fprintf('----------------------------------\n');
fprintf('| Processing Stage   | SNR (dB)   |\n');
fprintf('|--------------------|------------|\n');
fprintf('| Original Noisy     | %7.2f     |\n', snr_original);
fprintf('| After Notch        | %7.2f     |\n', snr_notch);
fprintf('| After Wiener       | %7.2f     |\n', snr_wiener);
fprintf('| After LMS          | %7.2f     |\n', snr_lms);
fprintf('----------------------------------\n');

%% Audio Playback
% soundsc(noisy, fs); pause(length(noisy)/fs + 1)
% soundsc(lms_filtered, fs);

%% Helper Functions
function plot_fft(signal, fs, ax, title_str, color, linestyle)
    L = length(signal);
    f = fs*(0:(floor(L/2)))/L;
    Y = fft(signal);
    plot(ax, f, abs(Y(1:floor(L/2)+1)), 'Color', color, 'LineStyle', linestyle)
    title(ax, title_str), xlabel('Frequency (Hz)'), ylabel('Magnitude')
    xlim([0 fs/2]), grid(ax, 'on')
end

function [filtered, tonal_freqs] = adaptive_notch_filter(signal, fs)
    L = length(signal);
    Y = fft(signal);
    f = fs*(0:(floor(L/2)))/L;
    spectrum = abs(Y(1:floor(L/2)+1));
    
    [~,locs] = findpeaks(spectrum, 'MinPeakHeight',0.2*max(spectrum), 'MinPeakDistance',50);
    tonal_freqs = f(locs(f(locs) < 0.99*(fs/2)));
    
    filtered = signal;
    for k = 1:length(tonal_freqs)
        wo = tonal_freqs(k)/(fs/2);
        [b,a] = iirnotch(wo, 75/(fs/2));
        filtered = filter(b,a,filtered);
    end
end

function H = wiener_filter(noisy, noise)
    Pn = abs(fft(noise)).^2;
    Py = abs(fft(noisy)).^2;
    H = max(0, 1 - Pn./(Py + eps));
end

function lms_filtered = adaptive_lms(signal, noise_ref)
    order = 128;
    mu = 0.001;
    lms = dsp.LMSFilter(order, 'StepSize', mu, 'Method', 'Normalized LMS');
    [~, lms_filtered] = lms(noise_ref(1:length(signal)), signal);
end

function snr = calculate_snr(clean_sig, processed_sig)
    noise_component = clean_sig - processed_sig;
    signal_power = sum(clean_sig.^2);
    noise_power = sum(noise_component.^2);
    snr = 10 * log10(signal_power / noise_power);
end
