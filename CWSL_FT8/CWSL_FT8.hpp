
#include "unpack.h"
#include "ldpc.h"
#include "decode.h"
#include "constants.h"
#include "encode.h"

#include "wave.h"
#include "debug.h"
#include "kiss_fftr.h"

#define LOG_LEVEL   LOG_INFO

#include "decode_audio_buffer.h"
#include "dsp.h"
#include "ring_buffer.h"


constexpr int kMax_candidates = 2500;
constexpr int kLDPC_iterations = 125;

constexpr int kMax_decoded_messages = 100;
constexpr int kMax_message_length = 40;


constexpr double fsk_dev = 6.25f;
constexpr int sample_rate = 6000;
constexpr int num_bins = static_cast<int>(sample_rate / (2 * fsk_dev));

constexpr int block_size = 2 * num_bins;

constexpr int num_samples = 15.5 * sample_rate;

constexpr int num_blocks = (num_samples - (block_size / 2) - block_size) / block_size;

constexpr int N = num_bins;

