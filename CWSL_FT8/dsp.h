#pragma once

#include <cstdint>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif // !M_PI

template <typename T>
T hann_i(const int i, const int N) {
	const float N_f = static_cast<T>(N);
	const float i_f = static_cast<T>(i);
	const float x = std::sin(M_PI * i_f / (N_f - 1.0f));
	return x * x;
}

template <typename T>
T hamming_i(const int i, const int N) {
	const float N_f = static_cast<T>(N);
	const float i_f = static_cast<T>(i);

	constexpr float a0 = 25.0f / 46.0f;
	constexpr float a1 = 1.0f - a0;

	float x1 = std::cos(2.0f * M_PI * i_f / (N_f - 1));
	return a0 - a1 * x1;
}

template <typename T>
T blackman_i(const int i, const int N) {
	const float N_f = static_cast<T>(N);
	const float i_f = static_cast<T>(i);
	constexpr float alpha = 0.16f; // or 2860/18608
	constexpr float a0 = (1.0f - alpha) / 2.0f;
	constexpr float a1 = 1.0f / 2;
	constexpr float a2 = alpha / 2;

	const float x1 = std::cos(2.0f * M_PI * i_f / (N_f - 1.0f));
	const float x2 = 2 * x1 * x1 - 1; // Use double angle formula

	return a0 - a1 * x1 + a2 * x2;
}

template <typename T>
void normalize_signal(T* signal, const size_t num_samples) {
	constexpr T max_amp = 1E-5f;
	for (int i = 0; i < num_samples; ++i) {
		const T amp = std::abs(signal[i]);
		if (amp > max_amp) {
			max_amp = amp;
		}
	}
	for (int i = 0; i < num_samples; ++i) {
		signal[i] /= max_amp;
	}
}


// Compute FFT magnitudes (log power) for each timeslot in the signal
template <int B>
void extract_power(const float signal[], int num_blocks, uint8_t power[]) {
	constexpr int block_size = 2 * B;      // Average over 2 bins per FSK tone
	constexpr int nfft = 2 * block_size;          // We take FFT of two blocks, advancing by one
	constexpr float fft_norm = 2.0f / nfft;

	float   window[nfft];

	for (int i = 0; i < nfft; ++i) {
		window[i] = blackman_i<float>(i, nfft);
	}

	size_t  fft_work_size;
	kiss_fftr_alloc(nfft, 0, 0, &fft_work_size);

	LOG(LOG_INFO, "N_FFT = %d\n", nfft);
	LOG(LOG_INFO, "FFT work area = %lu\n", fft_work_size);

	void* fft_work = malloc(fft_work_size);
	kiss_fftr_cfg fft_cfg = kiss_fftr_alloc(nfft, 0, fft_work, &fft_work_size);

	int offset = 0;
	float max_mag = -100.0f;
	for (int i = 0; i < num_blocks; ++i) {
		// Loop over two possible time offsets (0 and block_size/2)
		for (int time_sub = 0; time_sub <= block_size / 2; time_sub += block_size / 2) {
			kiss_fft_scalar timedata[nfft];
			kiss_fft_cpx    freqdata[nfft / 2 + 1];
			float           mag_db[nfft / 2 + 1];

			// Extract windowed signal block
			for (int j = 0; j < nfft; ++j) {
				timedata[j] = window[j] * signal[(i * block_size) + (j + time_sub)];
			}

			kiss_fftr(fft_cfg, timedata, freqdata);

			// Compute log magnitude in decibels
			for (int j = 0; j < nfft / 2 + 1; ++j) {
				float mag2 = (freqdata[j].i * freqdata[j].i + freqdata[j].r * freqdata[j].r);
				mag_db[j] = 10.0f * log10f(1E-10f + mag2 * fft_norm * fft_norm);
			}

			// Loop over two possible frequency bin offsets (for averaging)
			for (int freq_sub = 0; freq_sub < 2; ++freq_sub) {
				for (int j = 0; j < B; ++j) {
					float db1 = mag_db[j * 2 + freq_sub];
					float db2 = mag_db[j * 2 + freq_sub + 1];
					float db = (db1 + db2) / 2;

					// Scale decibels to unsigned 8-bit range and clamp the value
					int scaled = (int)(2 * (db + 120));
					power[offset] = (scaled < 0) ? 0 : ((scaled > 255) ? 255 : scaled);
					++offset;

					if (db > max_mag) max_mag = db;
				}
			}
		}
	}

	LOG(LOG_INFO, "Max magnitude: %.1f dB\n", max_mag);
	free(fft_work);
	free(window);
}
