import numpy as np
from scipy.signal import tukey

# === Parameters ===
# input_file = 'noisy_speech.txt'
input_file = 'external_noise.txt'
# output_file = 'noisy_speech_toned.txt'
output_file = 'external_noise_toned.txt'
sampling_rate = 44100  # Adjust if needed

# === Load original signal ===
signal1 = np.loadtxt(input_file)
signal = np.zeros_like(signal1)

# === Tone definitions: (start_time, end_time, frequency in Hz, magnitude in dB) ===
tones = [
    (3, 5, 2000, 0.4),
    # (4, 8, 5000, 0.2),
    # (8, 9, 300, 0.5),
    # (5, 9, 10000, 0.1)
]

# === Add tones with smoothing ===
for start_time, end_time, freq, db in tones:
    # Sample range
    start_sample = int(start_time * sampling_rate)
    end_sample = int(end_time * sampling_rate)
    num_samples = end_sample - start_sample
    t = np.arange(num_samples) / sampling_rate

    # dB to linear scale (assuming 0 dB = full scale = 1.0)
    magnitude = 10 ** (db / 20)

    # Generate the sine wave
    tone = magnitude * np.sin(2 * np.pi * freq * t)

    # Apply a Tukey window to smooth start/end and reduce spectral leakage
    # window = tukey(num_samples, alpha=0.5)  # alpha=0.5 = moderate smoothing
    # tone *= window

    # Add tone into the original signal
    signal[start_sample:end_sample] += tone

# === Save the modified signal ===
np.savetxt(output_file, signal)
