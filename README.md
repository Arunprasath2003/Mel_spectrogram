# Mel_spectrogram
<br>

**FFT Implementation:**
<br>

* Utilizes the Cooley-Tukey algorithm to perform Fast Fourier Transform (FFT).
* The fft function is responsible for recursively computing the FFT on a given vector of complex numbers.
<br>

**STFT Computation:**

<br>

* The stft function calculates the STFT of an audio signal.
* It uses a Hamming window function to apply windowing to the audio signal before computing the FFT.
* The result is a 2D vector representing the complex values of the STFT for different frames and frequencies.

<br>

**Mel Filterbank Calculation:**
<br>

* The mel_filter function applies Mel filterbank to the magnitude spectrum obtained from the STFT.
* Computes a triangular filterbank with a specified number of filters.
* The result is a 2D vector representing the Mel spectrogram.

<br>

Audio Loading:

<br>

* Utilizes the loadAudio function to load an audio file using the libsndfile library.
* Reads audio data and retrieves sample rate and channel information.

<br>

Main Function:

<br>

* Loads an audio file ("arun.wav").
* Defines parameters for STFT and Mel spectrogram computation (window size, hop size, and number of Mel filters).
* Computes the STFT and Mel spectrogram using the implemented functions.
* The magnitude spectrum is extracted from the STFT result for further analysis.
