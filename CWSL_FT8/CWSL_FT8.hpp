
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


constexpr int kMax_candidates = 100;
constexpr int kLDPC_iterations = 20;

constexpr int kMax_decoded_messages = 50;
constexpr int kMax_message_length = 20;


constexpr float fsk_dev = 6.25f;
constexpr int sample_rate = 48000;
constexpr int num_bins = static_cast<int>(sample_rate / (2 * fsk_dev));

const int N = num_bins;