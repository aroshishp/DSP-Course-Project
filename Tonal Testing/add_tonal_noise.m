function [noisy_speech_mod, external_noise_mod] = add_tonal_noise(noisy_speech, external_noise, freq, amp, fs)
    t = (0:length(noisy_speech)-1) / fs;
    tone = amp * sin(2 * pi * freq * t)';
    noisy_speech_mod = noisy_speech + tone;
    external_noise_mod = external_noise + tone;
end
