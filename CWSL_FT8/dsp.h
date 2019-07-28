#pragma once

#include <cstdint>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif // !M_PI

template <typename T>
T hann_i(const int i, const int N) {
	const double N_f = static_cast<T>(N);
	const double i_f = static_cast<T>(i);
	const double x = std::sin(M_PI * i_f / (N_f - 1.0f));
	return x * x;
}

template <typename T>
T hamming_i(const int i, const int N) {
	const double N_f = static_cast<T>(N);
	const double i_f = static_cast<T>(i);

	constexpr double a0 = 25.0f / 46.0f;
	constexpr double a1 = 1.0f - a0;

	double x1 = std::cos(2.0f * M_PI * i_f / (N_f - 1));
	return a0 - a1 * x1;
}

template <typename T>
T blackman_i(const int i, const int N) {
	const double N_f = static_cast<T>(N);
	const double i_f = static_cast<T>(i);
	constexpr double alpha = 0.16f; // or 2860/18608
	constexpr double a0 = (1.0f - alpha) / 2.0f;
	constexpr double a1 = 1.0f / 2;
	constexpr double a2 = alpha / 2;

	const double x1 = std::cos(2.0f * M_PI * i_f / (N_f - 1.0f));
	const double x2 = 2 * x1 * x1 - 1; // Use double angle formula

	return a0 - a1 * x1 + a2 * x2;
}

template <typename T>
T blackman_harris_i(const int i, const int N) {
    constexpr T a0 = 0.35875;
    constexpr T a1 = 0.48829;
    constexpr T a2 = 0.14128;
    constexpr T a3 = 0.01168;
    const T N_f = static_cast<T>(N);
    const T i_f = static_cast<T>(i);
    return a0 - (a1 * std::cos((2.0f * M_PI * i_f) / (N_f - 1.0))) + (a2 * std::cos((4.0 * M_PI * i_f) / (N_f - 1.0))) - (a3 * std::cos((6.0 * M_PI * i_f) / (N_f - 1.0)));
}

template <typename T>
T cosine4_i(const int i, const int N) {
    constexpr T a0 = 3.635819267707608e-001;
    constexpr T a1 = 4.891774371450171e-001;
    constexpr T a2 = 1.365995139786921e-001;
    constexpr T a3 = 1.064112210553003e-002;
    const T i_f = static_cast<T>(i);
    const T N_f = static_cast<T>(N);
    const T w = a0 - a1 * std::cos(2.0 * M_PI * i_f / (N_f - 1.0)) + a2 * cos(4.0 * M_PI * i_f / (N_f - 1.0)) - a3 * cos(6.0 * M_PI * i_f / (N_f - 1.0));
    return w;
}

template <typename T>
T cosine11_i(const int i, const int N) {
    constexpr T a0 = 2.151527506679809e-001;
    constexpr T a1 = 3.731348357785249e-001;
    constexpr T a2 = 2.424243358446660e-001;
    constexpr T a3 = 1.166907592689211e-001;
    constexpr T a4 = 4.077422105878731e-002;
    constexpr T a5 = 1.000904500852923e-002;
    constexpr T a6 = 1.639806917362033e-003;
    constexpr T a7 = 1.651660820997142e-004;
    constexpr T a8 = 8.884663168541479e-006;
    constexpr T a9 = 1.938617116029048e-007;
    constexpr T a10 = 8.482485599330470e-010;
    const T i_f = static_cast<T>(i);
    const T N_f = static_cast<T>(N);

    const T w = a0 - a1 * std::cos(2.0 * M_PI * i_f / (N_f - 1.0)) + a2 * std::cos(4.0 * M_PI * i_f / (N_f - 1.0)) - a3 * std::cos(6.0 * M_PI * i_f / (N_f - 1.0)) + a4 * std::cos(8.0 * M_PI * i_f / (N_f - 1.0)) - a5 * std::cos(10.0 * M_PI * i_f / (N_f - 1.0)) + a6 * std::cos(12.0 * M_PI * i_f / (N_f - 1.0)) - a7 * std::cos(14.0 * M_PI * i_f / (N_f - 1.0)) + a8 * std::cos(16.0 * M_PI * i_f / (N_f - 1.0)) - a9 * std::cos(18.0 * M_PI * i_f / (N_f - 1.0)) + a10 * std::cos(20.0 * M_PI * i_f / (N_f - 1.0));
    return w;
}

template <typename T,int N>
void normalize_signal(std::array<T,N>& signal) {
	const T max_amp = *std::max_element(signal.begin(),signal.end());
	for (size_t i = 0; i < signal.size(); ++i) {
		signal[i] /= max_amp;
	}
}

// Compute FFT magnitudes (log power) for each timeslot in the signal
template <int B, int num_samples>
void extract_power(const std::array<double, num_samples>& signal, int num_blocks, std::vector<uint8_t>& power) {
	constexpr int block_size = 2 * B;      // Average over 2 bins per FSK tone
	constexpr int nfft = 2 * block_size;          // We take FFT of two blocks, advancing by one
	constexpr double fft_norm = 2.0f / nfft;

	std::array<double, nfft> window1;

	for (int i = 0; i < nfft; ++i) {
		window1[i] = blackman_harris_i<double>(i, nfft);
	}

    std::array<double, nfft / 2 + 1> window2;

    for (int i = 0; i < nfft / 2 + 1; ++i) {
        window2[i] = cosine11_i<double>(i, nfft / 2 + 1);
    }

	size_t  fft_work_size;
	kiss_fftr_alloc(nfft, 0, 0, &fft_work_size);

	//LOG(LOG_INFO, "N_FFT = %d\n", nfft);
	//LOG(LOG_INFO, "FFT work area = %lu\n", fft_work_size);

	void* fft_work = malloc(fft_work_size);
	kiss_fftr_cfg fft_cfg = kiss_fftr_alloc(nfft, 0, fft_work, &fft_work_size);

	int offset = 0;
	double max_mag = -INFINITY;
	for (int i = 0; i < num_blocks; ++i) {
		// Loop over two possible time offsets (0 and block_size/2)
		for (int time_sub = 0; time_sub <= block_size / 2; time_sub += block_size / 2) {
			kiss_fft_scalar timedata[nfft];
			kiss_fft_cpx    freqdata[nfft / 2 + 1];
			double           mag_db[nfft / 2 + 1];

			// Extract windowed signal block
			for (int j = 0; j < nfft; ++j) {
				timedata[j] = window1[j] * signal[(i * block_size) + (j + time_sub)];
                //timedata[j] = signal[(i * block_size) + (j + time_sub)];

			}

			kiss_fftr(fft_cfg, timedata, freqdata);

			// Compute log magnitude in decibels
			for (int j = 0; j < nfft / 2 + 1; ++j) {
				double mag2 = (freqdata[j].i * freqdata[j].i + freqdata[j].r * freqdata[j].r);

            //    double mag2 = (window2[j]*freqdata[j].i * window2[j]*freqdata[j].i + window2[j]*freqdata[j].r * window2[j]*freqdata[j].r);


				mag_db[j] = 10.0 * std::log10(1E-10 + mag2 * fft_norm * fft_norm);
			}

			// Loop over two possible frequency bin offsets (for averaging)
			for (int freq_sub = 0; freq_sub < 2; ++freq_sub) {
				for (int j = 0; j < B; ++j) {
					double db1 = mag_db[j * 2 + freq_sub];
					double db2 = mag_db[j * 2 + freq_sub + 1];
					double db = (db1 + db2) / 2;

					// Scale decibels to unsigned 8-bit range and clamp the value
					int scaled = (int)(2 * (db + 120));
					power[offset] = (scaled < 0) ? 0 : ((scaled > 255) ? 255 : scaled);
					++offset;

					if (db > max_mag) max_mag = db;
				}
			}
		}
	}

	//LOG(LOG_INFO, "Max magnitude: %.1f dB\n", max_mag);
	free(fft_work);
}
