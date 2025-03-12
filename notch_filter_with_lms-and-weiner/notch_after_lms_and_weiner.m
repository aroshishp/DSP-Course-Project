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
figure('Name','Notch After Processing', 'Position', [100 100 1400 1500])

%% Stage 0: Clean Signal Reference
subplot(6,2,1)
plot_fft(clean, fs, gca, 'Clean Reference', 'b', '-')
subplot(6,2,2)
plot(t, clean)
title('Time Domain - Clean Speech')
xlabel('Time (s)'), axis tight

%% Stage 1: Noisy Signal
subplot(6,2,3)
plot_fft(noisy, fs, gca, 'Original Noisy', 'r', '-')
subplot(6,2,4)
plot(t, noisy)
title('Time Domain - Noisy Speech')
xlabel('Time (s)'), axis tight

%% Stage 2: Wiener Filtering
H_wiener = (abs(fft(clean)).^2) ./ (abs(fft(clean)).^2 + abs(fft(noise)).^2);
wiener_filtered = real(ifft(H_wiener .* fft(noisy)));

subplot(6,2,5)
plot_fft(wiener_filtered, fs, gca, 'After Wiener', 'g', '-')
hold on
plot_fft(clean, fs, gca, '', 'b', ':')
hold off
subplot(6,2,6)
plot(t, wiener_filtered)
title('Time Domain - Wiener Filtered')
xlabel('Time (s)'), axis tight

%% Stage 3: LMS Filtering
order = 128;
mu = 0.001;
lms = dsp.LMSFilter(order, 'StepSize', mu, 'Method', 'Normalized LMS');
[~, lms_filtered] = lms(noise(1:length(wiener_filtered)), wiener_filtered);

subplot(6,2,7)
plot_fft(lms_filtered, fs, gca, 'After LMS', 'm', '-')
hold on
plot_fft(clean, fs, gca, '', 'b', ':')
hold off
subplot(6,2,8)
plot(t, lms_filtered)
title('Time Domain - LMS Filtered')
xlabel('Time (s)'), axis tight

%% Stage 4: Notch on Wiener Output
[notch_wiener, ~] = apply_notch(wiener_filtered, fs);

subplot(6,2,9)
plot_fft(notch_wiener, fs, gca, 'Wiener + Notch', 'c', '-')
hold on
plot_fft(clean, fs, gca, '', 'b', ':')
hold off
subplot(6,2,10)
plot(t, notch_wiener)
title('Time Domain - Wiener Notched')
xlabel('Time (s)'), axis tight

%% Stage 5: Notch on LMS Output
[notch_lms, ~] = apply_notch(lms_filtered, fs);

subplot(6,2,11)
plot_fft(notch_lms, fs, gca, 'LMS + Notch', 'k', '-')
hold on
plot_fft(clean, fs, gca, '', 'b', ':')
hold off
subplot(6,2,12)
plot(t, notch_lms)
title('Time Domain - LMS Notched')
xlabel('Time (s)'), axis tight

%% SNR Calculations
snr_results = zeros(6,1);
snr_results(1) = calculate_snr(clean, noisy);
snr_results(2) = calculate_snr(clean, wiener_filtered);
snr_results(3) = calculate_snr(clean, lms_filtered);
snr_results(4) = calculate_snr(clean, notch_wiener);
snr_results(5) = calculate_snr(clean, notch_lms);

fprintf('\nSNR Improvement Report:\n');
fprintf('-----------------------------------------\n');
fprintf('| Processing Stage    | SNR (dB) |\n');
fprintf('|---------------------|----------|\n');
fprintf('| Original Noisy      | %7.2f  |\n', snr_results(1));
fprintf('| After Wiener        | %7.2f  |\n', snr_results(2));
fprintf('| After LMS           | %7.2f  |\n', snr_results(3));
fprintf('| Wiener + Notch      | %7.2f  |\n', snr_results(4));
fprintf('| LMS + Notch         | %7.2f  |\n', snr_results(5));
fprintf('-----------------------------------------\n');

%% Audio Playback
soundsc(noisy, fs); pause(length(noisy)/fs + 1)
soundsc(wiener_filtered, fs); pause(length(wiener_filtered)/fs + 1)
soundsc(lms_filtered, fs); pause(length(lms_filtered)/fs + 1)
soundsc(notch_wiener, fs); pause(length(notch_wiener)/fs + 1)
soundsc(notch_lms, fs);

%% Notch Filter Function
function [filtered, tonal_freqs] = apply_notch(input_sig, fs)
    Y = fft(input_sig);
    L = length(input_sig);
    f = fs*(0:(floor(L/2)))/L;
    spectrum = abs(Y(1:floor(L/2)+1));
    
    threshold = 0.2 * max(spectrum);
    [~, locs] = findpeaks(spectrum, 'MinPeakHeight', threshold, 'MinPeakDistance', 50);
    tonal_freqs = f(locs);
    valid_freqs = tonal_freqs(tonal_freqs < fs/2*0.99);

    filtered = input_sig;
    for k = 1:length(valid_freqs)
        wo = valid_freqs(k)/(fs/2);
        bw = 75/(fs/2);
        [b,a] = iirnotch(wo, bw);
        filtered = filter(b,a,filtered);
    end
end

%% FFT Plotting Function
function plot_fft(signal, fs, ax, title_str, color, linestyle)
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
