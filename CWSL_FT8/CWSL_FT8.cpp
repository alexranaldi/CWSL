// Alex Ranaldi  W2AXR   alexranaldi@gmail.com

// LICENSE: GNU General Public License v3.0
// THE SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.

// This work is based on: https://github.com/HrochL/CWSL

// This program was written to feed audio from QS1R SDR to WSJT-X for FT8.
//  Virtual Audio Cables are used to pass audio.

#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <string>
#include <cstdint>
#include <vector>
#include <thread>
#include <chrono>
#include <atomic>
#include <ctime>
#include <iomanip>

#include "../Utils/SharedMemory.h"

// SSB demod
#include "../Utils/SSBD.hpp"
// Auto scale audio
#include "../Utils/AutoScaleAF.hpp"
// Up sample audio
#include "../Utils/Upsampler.hpp"

#include "CWSL_FT8.hpp"


// Maximum of CWSL bands
#define MAX_CWSL   32

const int BITS_PER_SAMPLE = 32;

// Prefix and suffix for the shared memories names
const std::string gPreSM = "CWSL";
std::string gPostSM = "Band";

std::atomic_bool terminateFlag;
std::atomic_bool holdScaleFactor;


ring_buffer_t<std::complex<float>*> iq_buffer;
decode_audio_buffer_t<float, num_samples> decode_audio_buffer;


int SF = 16;

int SMNumber = -1;

std::string createSharedMemName(const int bandIndex) {
	// create name of shared memory
	std::string Name = gPreSM + std::to_string(bandIndex) + gPostSM;
	if (SMNumber != -1) {
		Name += std::to_string(SMNumber);
	}
	return Name;
}

void waitForTime() {
	while (!terminateFlag) {
		std::time_t t = std::time(nullptr);
		const auto ts = std::gmtime(&t);
	//	std::cout << "UTC:   " << std::put_time(ts, "%c %Z") << std::endl;
		const bool go = ts->tm_sec == 0 || ts->tm_sec == 15 || ts->tm_sec == 30 || ts->tm_sec == 45;
		if (go) {
			return;
		}
		else {
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		}
	}
}
	
void doDecodeFT8() {
	std::cout << "Decoding..." << std::endl;

	// Compute FFT over the whole signal and store it
	uint8_t power[num_blocks * 4 * num_bins];
	extract_power<num_bins>(decode_audio_buffer.buf.data(), num_blocks, power);
	
	// Find top candidates by Costas sync score and localize them in time and frequency
	ft8::Candidate candidate_list[kMax_candidates];
	int num_candidates = ft8::find_sync(power, num_blocks, num_bins, ft8::kCostas_map, kMax_candidates, candidate_list);

	// Go over candidates and attempt to decode messages
	char    decoded[kMax_decoded_messages][kMax_message_length];
	int     num_decoded = 0;
	for (int idx = 0; idx < num_candidates; ++idx) {
		ft8::Candidate& cand = candidate_list[idx];
		float freq_hz = (cand.freq_offset + cand.freq_sub / 2.0f) * fsk_dev;
		float time_sec = (cand.time_offset + cand.time_sub / 2.0f) / fsk_dev;

		float   log174[ft8::N];
		ft8::extract_likelihood(power, num_bins, cand, ft8::kGray_map, log174);

		// bp_decode() produces better decodes, uses way less memory
		uint8_t plain[ft8::N];
		int     n_errors = 0;
		ft8::bp_decode(log174, kLDPC_iterations, plain, &n_errors);
		//ldpc_decode(log174, kLDPC_iterations, plain, &n_errors);

		if (n_errors > 0) {
			LOG(LOG_DEBUG, "ldpc_decode() = %d (%.0f Hz)\n", n_errors, freq_hz);
			continue;
		}

		// Extract payload + CRC (first ft8::K bits)
		uint8_t a91[ft8::K_BYTES];
		ft8::pack_bits(plain, ft8::K, a91);

		// Extract CRC and check it
		uint16_t chksum = ((a91[9] & 0x07) << 11) | (a91[10] << 3) | (a91[11] >> 5);
		a91[9] &= 0xF8;
		a91[10] = 0;
		a91[11] = 0;
		uint16_t chksum2 = ft8::crc(a91, 96 - 14);
		if (chksum != chksum2) {
			LOG(LOG_DEBUG, "Checksum: message = %04x, CRC = %04x\n", chksum, chksum2);
			continue;
		}

		char message[kMax_message_length];
		ft8::unpack77(a91, message);

		// Check for duplicate messages (TODO: use hashing)
		bool found = false;
		for (int i = 0; i < num_decoded; ++i) {
			if (0 == strcmp(decoded[i], message)) {
				found = true;
				break;
			}
		}

		if (!found && num_decoded < kMax_decoded_messages) {
			strcpy(decoded[num_decoded], message);
			++num_decoded;

			// Fake WSJT-X-like output for now
			int snr = 0;    // TODO: compute SNR
			printf("000000 %3d %4.1f %4d ~  %s\n", cand.score, time_sec, (int)(freq_hz + 0.5f), message);
		}
	}
	LOG(LOG_INFO, "Decoded %d messages\n", num_decoded);



	std::cout << "Done decoding" << std::endl;
}

void decodeLoop() {
	while (!terminateFlag) {
		waitForTime();
		std::time_t t = std::time(nullptr);
		const auto ts = std::gmtime(&t);
		std::cout << "UTC:   " << std::put_time(ts, "%c %Z") << std::endl;
		if (terminateFlag) {
			return;
		}
		doDecodeFT8();
		if (terminateFlag) {
			return;
		}
		decode_audio_buffer.clear();
		// wait at least 1.5 seconds so we don't double decode
		std::this_thread::sleep_for(std::chrono::milliseconds(1500));
	}
}

void readIQ(CSharedMemory &SM, const size_t iq_len) {

    while (!terminateFlag) {
    
        // wait for new data from receiver. Blocks until data received
        SM.WaitForNewData();
        
        // wait for a slot to be ready in the buffer. Blocks until slot available.
        if (!iq_buffer.wait_for_empty_slot()) {
            std::cout << "No slots available in IQ buffer!" << std::endl;
            continue;
        }

        std::complex<float>* iq_data_f = iq_buffer.recs[iq_buffer.write_index];

        // read block of data from receiver
        const bool readSuccess = SM.Read((PBYTE)iq_data_f, iq_len * sizeof(std::complex<float>));
        if (readSuccess) {
			iq_buffer.inc_write_index();
        }
        else {
            std::cout << "Did not read any I/Q data from shared memory" << std::endl;
        }
    }
}

template <typename T, int N>
void demodulate(SSBD<T>& ssbd, Upsampler<T>& upsamp, AutoScaleAF<T>& af, const size_t iq_len, const size_t decRatio, decode_audio_buffer_t<T,N>& decode_audio_buffer) {
    // Demodulate IQ for SSB

    const size_t ssbd_in_size = ssbd.GetInSize();
    const size_t upsamp_ratio = upsamp.GetRatio();

    std::vector<T> af6khz(iq_len / decRatio, 0.0);

    const size_t n48khz = iq_len * upsamp_ratio / decRatio;

    T scaleFactor = static_cast<T>(std::pow(SF, 2));
    T lastScaleFactor = 0;

	std::vector<T> af48khz(n48khz);

    while (!terminateFlag) {
        // get IQ data
        std::complex<T> *xc = iq_buffer.pop();

        for (size_t n = 0; n < iq_len; n += ssbd_in_size) {
            ssbd.Iterate(xc + n, af6khz.data() + n / decRatio);
        }

        // Upsample 6kHz audio to 48khz
        for (size_t n = 0; n < af6khz.size(); ++n) {
            upsamp.Iterate(af6khz.data() + n, af48khz.data() + n * upsamp_ratio);
        }

        // If auto AF is enabled
        if (-1 == SF && !holdScaleFactor) {
            scaleFactor = af.getScaleFactor(af48khz);
            if (scaleFactor != lastScaleFactor) {
                std::cout << "Scale factor set to: " << scaleFactor << std::endl;
                lastScaleFactor = scaleFactor;
            }
        }
        for (size_t n = 0; n < n48khz; ++n) {
            af48khz[n] *= scaleFactor;
        }

        // AF write complete
		decode_audio_buffer.write(af48khz);
    }
}


//////////////////////////////////////////////////////////////////////////////
// Find the right band
int findBand(const int64_t f) {
    CSharedMemory SM;
    SM_HDR h;

    // try to find right band - for all possible bands ...
    for (int bandIndex = 0; bandIndex < MAX_CWSL; ++bandIndex) {

        // create name of shared memory
		const std::string Name = createSharedMemName(bandIndex);

        // try to open shared memory
        if (SM.Open(Name.c_str())) {
            // save data from header of this band
            memcpy(&h, SM.GetHeader(), sizeof(SM_HDR));

            // close shared memory
            SM.Close();

            // is frequency in this band ?
            if ((h.SampleRate > 0) && (f >= h.L0 - h.SampleRate / 2) && (f <= h.L0 + h.SampleRate / 2)) {
                // yes -> assign it and break the finding loop
                return bandIndex;
            }
        }
    }
    return -1;
}


///////////////////////////////////////////////////////////////////////////////
// Main function
int main(int argc, char **argv)
{

    // check number of input parameters
    if (argc < 4) {
        // print usage
        std::cout << "Not enough input arguments!" << std::endl;
        std::cout << "Usage: CWSL_FT8 FreqHz Scale_factor Shared_Mem" << std::endl;
        std::cout << "    FreqHz is the frequency in Hz" << std::endl
                  << "    Scale_factor is -1 for Auto-Scaling" << std::endl
                  << "    Shared_Mem is an optional single numeric digit specifying the shared memory interface" << std::endl;
        std::cout << std::endl;           
        return EXIT_SUCCESS;
    }

    CSharedMemory SM;
    SM_HDR* SHDR;
    int nMem = 0;
    size_t nWO = 0;

    terminateFlag = false;

    // Shared Mem
    if (argc >= 4) {
        std::cout << "A shared memory interface was specified." << std::endl;
        if ((sscanf(argv[3], "%d", &SMNumber) != 1)) {
            std::cout << "Unable to deciper the specified Shared_Mem interface, so ignoring" << std::endl;
        }
        else {
            std::cout << "Using shared mem interface number: " << SMNumber << std::endl;
        }
    }

    int64_t ssbFreq = 0;

    // Get USB frequency
    if ((sscanf(argv[1], "%I64d", &ssbFreq) != 1)) {
        std::cout << "Unable to parse Frequency input" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "SSB Demodulator F=" << ssbFreq << std::endl;
    nMem = findBand(ssbFreq);
    if (-1 == nMem) {
        std::cout << "Bad FreqHz or Shared_Mem specified" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Using Receiver Band Index = " << nMem << std::endl;

    // Scale Factor
    if (argc >= 3) {
        if ((sscanf(argv[2], "%d", &SF) != 1) || (SF > 24)) {
            std::cout << "Bad Scale_factor: " << argv[2] << std::endl;
            return EXIT_FAILURE;
        }
    }
    else {
        // Use AutoAF
        SF = -1;
    }
    if (-1 == SF) {
        std::cout << "Using automatic scale factor (AutoAF)" << std::endl;
    }
    else {
        std::cout << "Using user-specified scale factor: " << SF << std::endl;
    }

    // USB/LSB.  USB = 1, LSB = 0
    const int USB = 1;
 
    //
    // Setup shared memory and receiver interface
    //
   
    // try to open shared memory
	const std::string name = createSharedMemName(nMem);
    if (!SM.Open(name.c_str())) {
        fprintf(stderr, "Can't open shared memory for %d receiver\n", nMem);
        return EXIT_FAILURE;
    }
    // get info about channel
    SHDR = SM.GetHeader();
    const size_t SR = static_cast<size_t>(SHDR->SampleRate);
    const int BIS = SHDR->BlockInSamples;
    std::cout << "Receiver: " << nMem
        << "\tSample Rate: " << SR
        << "\tBlock In Samples: " << BIS
        << "\tLO: " << SHDR->L0
        << std::endl;

    //
    // Create the SSB demodulator
    //
  
    const size_t SSB_BW = 3000;
    std::cout << "SSB Bandwidth: " << SSB_BW << " Hz" << std::endl;
    // F is always Fc-LO
    const int64_t LO = static_cast<int64_t>(SHDR->L0);
    const int64_t F = ssbFreq - SHDR->L0;
    SSBD<float> ssbd(SR, SSB_BW, static_cast<float>(F), static_cast<bool>(USB));
    const size_t SSB_SR = ssbd.GetOutRate();
    const size_t decRatio = SR / SSB_SR;

    //
    // Create the UpSampler
    //

    Upsampler<float> upsamp(3); //2^3=8 ratio
    const size_t Wave_SR = 48000;
    std::cout << "Upsampling audio from " << SSB_SR << " Hz to " << Wave_SR << " Hz" << std::endl;
    const size_t upsamp_ratio = upsamp.GetRatio();
    if (upsamp_ratio != Wave_SR / SSB_SR) {
        std::cout << "UpSampler does not support the specified wave output sample rate" << std::endl;
        SM.Close();
        return EXIT_FAILURE;
    }

    // 
    // Open WinWave (audio device)
    //

    const size_t iq_len = BIS;
    const size_t af_len = iq_len * upsamp_ratio / decRatio;

    //
    // Create AutoAF. Only used if SF == -1
    //

    const float headroomdB = 18;
    const float windowdB = 22;
    const size_t numSampBeforeAfInc = Wave_SR * 300; // 300s of samples
	const float clipVal = static_cast<float>(std::pow(2, BITS_PER_SAMPLE - 1) - 1);

    AutoScaleAF<float> af(headroomdB, windowdB, clipVal, numSampBeforeAfInc);
    if (-1 == SF) {
        float autoAfhi = 0;
        float autoAflo = 0;
        af.getThresholds(autoAfhi, autoAflo);
        std::cout << "AutoAF ideal headroom [max,min] from clip: [" << autoAfhi << " dB, " << autoAflo << " dB]" << std::endl;
        std::cout << std::endl << "Press H to enable/disable hold of scale factor" << std::endl;
    }
    holdScaleFactor = false;

    //
    // Prepare circular buffers
    // 

    bool initStatus = iq_buffer.initialize(4096);
	if (!initStatus) {
		std::cout << "Failed to initialize memory" << std::endl;
		return EXIT_FAILURE;
	}
    for (size_t k = 0; k < iq_buffer.size; ++k) {
        iq_buffer.recs[k] = reinterpret_cast<std::complex<float>*>( malloc( sizeof(std::complex<float>) * iq_len ) );
    }

    //
    //  Start Threads
    //

    std::cout << "Creating receiver thread..." << std::endl;
    std::thread iqThread = std::thread(readIQ, std::ref(SM), BIS);
    std::cout << "Creating SSB Demodulator thread..." << std::endl;
    std::thread demodThread = std::thread(&demodulate<float,num_samples>, std::ref(ssbd), std::ref(upsamp), std::ref(af), BIS, decRatio, std::ref(decode_audio_buffer));
	std::cout << "Creating decoder thread..." << std::endl;
	std::thread decodeThread = std::thread(&decodeLoop);

	/*
    std::cout << "Creating Audio Device thread..." << std::endl;
    std::thread audioThread = std::thread(&writeAudio, std::ref(wave), af_len);
    */

    //
    //  Main Loop
    //

    std::cout << std::endl << "Main loop started! Press Q to terminate." << std::endl;
    while (!terminateFlag) {
        // Was Exit requested?
        if (_kbhit()) {
            const char ch = _getch();
            if (ch == 'Q' || ch == 'q') {
                std::cout << "Q pressed, so terminating" << std::endl;
                terminateFlag = true;
            }
            else if (SF == -1 && (ch == 'H' || ch == 'h')) {
                holdScaleFactor = !holdScaleFactor;
                if (holdScaleFactor){
                    std::cout << "Holding scale factor" << std::endl;
                }
                else{
                    std::cout << "No longer holding scale factor" << std::endl;
                }
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }

    //
    // Join all threads
    //
    
	/*
    if (audioThread.joinable()) {
        audioThread.join();
    }
	*/
    if (demodThread.joinable()) {
        demodThread.join();
    }
    if (iqThread.joinable()) {
        iqThread.join();
    }
    


    //
    //  Clean up
    //

    //wave.stop();
    SM.Close();

    std::cout << "Deallocating memory..." << std::endl;

    // free memory
  

    return EXIT_SUCCESS;

}    