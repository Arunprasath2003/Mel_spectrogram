import torch
import torchaudio
from stft_mel_py import stft_mel_python
from stft_mel_cpp import stft, mel_filter

def test_stft_mel():
    # Load audio using torchaudio
    waveform, sample_rate = torchaudio.load('arun.wav', normalize=True)

    # Parameters
    window_size = 400
    hop_size = 160
    num_mel_filters = 80

    # Python implementation
    magnitude_py, _, mel_spectrogram_py = stft_mel_python(waveform, window_size, hop_size, num_mel_filters)

    # C++ implementation
    magnitude_cpp = stft(waveform.numpy().tolist(), window_size, hop_size)
    mel_spectrogram_cpp = mel_filter(magnitude_cpp, num_mel_filters)

    # Compare the results
    assert torch.allclose(torch.tensor(magnitude_cpp), magnitude_py, atol=1e-5)
    assert torch.allclose(torch.tensor(mel_spectrogram_cpp), mel_spectrogram_py, atol=1e-5)


