#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <sndfile.h>

// FFT Class
class FFT {
public:
    // Function to perform FFT (Cooley-Tukey algorithm)
    static void compute(std::vector<std::complex<double>>& a, bool inverse = false) {
        int n = a.size();
        if (n <= 1) return;

        // Divide
        std::vector<std::complex<double>> a0(n / 2), a1(n / 2);
        for (int i = 0; i < n / 2; i++) {
            a0[i] = a[i * 2];
            a1[i] = a[i * 2 + 1];
        }

        // Recursive FFT
        compute(a0, inverse);
        compute(a1, inverse);

        // Combine
        double theta = (inverse ? 1 : -1) * 2.0 * M_PI / n;
        /*
        std::complex<double> w(1): This initializes a complex number w with a real part of 1 and an imaginary part of 0. In the context of FFT (Fast Fourier 			Transform), w is commonly used as the twiddle factor, representing the rotation in the complex plane during the FFT computation.

		wn(cos(theta), sin(theta));: This part initializes another complex number wn using the polar form. The real part is cos(theta) and the imaginary part is 			sin(theta). The values of cos(theta) and sin(theta) are calculated based on the angle theta which is determined by the FFT direction (forward or inverse) 			and the current stage of the FFT algorithm.

		So, together, these two complex numbers w and wn are used in the FFT algorithm to perform the rotation of complex values during the divide-and-conquer 			process, contributing to the computation of the FFT.
        */
        std::complex<double> w(1), wn(cos(theta), sin(theta));
		//The following loop combines the results of the recursive FFTs:
		/*
		It iterates over the first half of the input sequence (a).
		t is the complex exponential factor multiplied by the corresponding element from the second half (a1).
		The result is a combination of elements from the first half (a0) and the rotated second half (a1).
		The loop updates the values in the original sequence (a).
		*/
        for (int i = 0; i < n / 2; i++) {
    		// Multiply the twiddle factor 'w' with the corresponding element from the second half 'a1'
    		std::complex<double> t = w * a1[i];

    		// Butterfly operation: Combine the results from the first half 'a0' and the multiplied twiddle factor
    		a[i] = a0[i] + t;

    		// Butterfly operation: Subtract the results from the first half 'a0' and the multiplied twiddle factor
    		a[i + n / 2] = a0[i] - t;

    		// Update the twiddle factor for the next iteration
    		w *= wn;
		}

    }
};

// STFT Class
class STFT {
public:
    // Function to compute the STFT of an audio signal
    static std::vector<std::vector<std::complex<double>>> compute(const std::vector<double>& audio, int window_size, int hop_size) {
        int N = audio.size();
        int num_frames = (N - window_size) / hop_size + 1;

        // Hamming window function
        auto hamming = [](int n, int N) {
            return 0.54 - 0.46 * cos(2.0 * M_PI * n / (N - 1));
        };

        // Compute STFT
        std::vector<std::vector<std::complex<double>>> stft_result(num_frames, std::vector<std::complex<double>>(window_size));

        for (int i = 0; i < num_frames; ++i) {
            for (int j = 0; j < window_size; ++j) {
                int index = i * hop_size + j;
                stft_result[i][j] = audio[index] * hamming(j, window_size);
            }

            // Apply FFT
            FFT::compute(stft_result[i]);
        }

        return stft_result;
    }
};

// MelFilter Class
class MelFilter {
public:
    // Function to apply Mel-filter to the magnitude spectrum
    static std::vector<std::vector<double>> apply(const std::vector<std::vector<std::complex<double>>>& stft_result, int num_mel_filters) {
        // Mel filterbank
        std::vector<std::vector<double>> mel_filterbank(num_mel_filters, std::vector<double>(stft_result[0].size()));

        // Compute Mel filterbank
        for (int i = 0; i < num_mel_filters; ++i) {
            for (int j = 0; j < stft_result[0].size(); ++j) {
                // Dummy triangular filter shape
                mel_filterbank[i][j] = std::max(0.0, 1 - std::abs((j - (stft_result[0].size() - 1) / 2.0) / ((stft_result[0].size() - 1) / 2.0 - 1) - i / (num_mel_filters - 1.0)));
            }
        }

        // Apply Mel filterbank to magnitude spectrum
        std::vector<std::vector<double>> mel_spectrogram(stft_result.size(), std::vector<double>(num_mel_filters));

        for (int i = 0; i < stft_result.size(); ++i) {
            for (int j = 0; j < num_mel_filters; ++j) {
                double sum = 0.0;
                for (int k = 0; k < stft_result[0].size(); ++k) {
                    sum += std::abs(stft_result[i][k]) * mel_filterbank[j][k];
                }
                mel_spectrogram[i][j] = sum;
            }
        }

        return mel_spectrogram;
    }
};

// AudioLoader Class
class AudioLoader {
public:
    // Function to load audio from a WAV file
    static std::vector<double> load(const char* filename, int& sampleRate, int& numChannels) {
        SF_INFO sfinfo;
        SNDFILE* file = sf_open(filename, SFM_READ, &sfinfo);
        if (!file) {
            std::cerr << "Error opening file: " << sf_strerror(NULL) << std::endl;
            return {};
        }

        // Get audio file information
        numChannels = sfinfo.channels;
        sampleRate = sfinfo.samplerate;

        // Read audio data
        std::vector<double> audio(sfinfo.frames * sfinfo.channels);
        sf_readf_double(file, audio.data(), sfinfo.frames);

        // Close the file
        sf_close(file);

        return audio;
    }
};

int main() {
    const char* filename = "arun.wav";
    int sampleRate, numChannels;

    // Load audio
    std::vector<double> audio = AudioLoader::load(filename, sampleRate, numChannels);

    // Parameters
    int window_size = 400;
    int hop_size = 160;
    int num_mel_filters = 80;

    // Compute STFT
    std::vector<std::vector<std::complex<double>>> stft_result = STFT::compute(audio, window_size, hop_size);

    // Compute Magnitude Spectrum
    std::vector<std::vector<double>> magnitude(stft_result.size(), std::vector<double>(stft_result[0].size()));
    for (int i = 0; i < stft_result.size(); ++i) {
        for (int j = 0; j < stft_result[i].size(); ++j) {
            magnitude[i][j] = std::abs(stft_result[i][j]);
        }
    }

    // Apply Mel filter
    std::vector<std::vector<double>> mel_spectrogram = MelFilter::apply(stft_result, num_mel_filters);

    return 0;
}

