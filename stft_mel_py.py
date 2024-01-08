import torch
import torchaudio
from torch.utils.cpp_extension import load

# Load C++ module using Pybind
stft_mel_cpp = load(name='stft_mel_cpp', sources=['stft_mel_cpp.cpp'])

def stft_mel_python(audio, window_size, hop_size, num_mel_filters):
    # STFT using PyTorch
    stft_result = torchaudio.transforms.Spectrogram()(audio)

    # Convert STFT to magnitude and phase
    magnitude = stft_result.abs()
    phase = torch.angle(stft_result)

    # Mel-filter using PyTorch
    mel_spectrogram = torchaudio.transforms.MelScale(num_mel_filters)(magnitude)

    return magnitude, phase, mel_spectrogram

# Load audio using torchaudio
waveform, sample_rate = torchaudio.load('arun.wav', normalize=True)

# Parameters
window_size = 400
hop_size = 160
num_mel_filters = 80

# Python implementation
magnitude_py, phase_py, mel_spectrogram_py = stft_mel_python(waveform, window_size, hop_size, num_mel_filters)

# C++ implementation
magnitude_cpp = stft_mel_cpp.stft(waveform.numpy().tolist(), window_size, hop_size)
mel_spectrogram_cpp = stft_mel_cpp.mel_filter(magnitude_cpp, num_mel_filters)

# Compare the results
assert torch.allclose(torch.tensor(magnitude_cpp), magnitude_py, atol=1e-5)
assert torch.allclose(torch.tensor(mel_spectrogram_cpp), mel_spectrogram_py, atol=1e-5)
print("Results match!")


