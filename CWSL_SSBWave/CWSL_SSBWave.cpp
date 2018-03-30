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

#include "../Utils/SharedMemory.h"

// SSB demod
#include "../Utils/SSBD.hpp"
// Auto scale audio
#include "../Utils/AutoScaleAF.hpp"
// Up sample audio
#include "../Utils/Upsampler.hpp"
// Write to Wave Out dev (virtual audio cable)
#include "../Utils/WinWave.hpp"

// Maximum of CWSL bands
#define MAX_CWSL   32

// Prefix and suffix for the shared memories names
const std::string gPreSM = "CWSL";
std::string gPostSM = "Band";

//////////////////////////////////////////////////////////////////////////////
// Find the right band
int findBand(const int64_t f) {
    CSharedMemory SM;
    SM_HDR h;
    int bandIndex;

    // try to find right band - for all possible bands ...
    for (bandIndex = 0; bandIndex < MAX_CWSL; ++bandIndex) {

        // create name of shared memory
        const std::string Name = gPreSM + std::to_string(bandIndex) + gPostSM;

        // try to open shared memory
        if (SM.Open(Name.c_str())) {
            // save data from header of this band
            memcpy(&h, SM.GetHeader(), sizeof(SM_HDR));

            // close shared memory
            SM.Close();

            // is frequeny into this band ?
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
    CSharedMemory SM;
    SM_HDR *SHDR;
    int SF = 16;
    int nMem = 0;
    size_t nWO = 0;

    WinWave wave = WinWave();

    // check number of input parameters
    if (argc < 4) {
        // print usage
        std::cout << "Not enough input arguments!" << std::endl;
        std::cout << "Usage: CWSL_SSBWave FreqHz WaveOutNr Scale_factor SSB_LSB" << std::endl;
        std::cout << "    FreqHz is the frequency in Hz" << std::endl
                  << "    WaveOutNr is the Wave Out device number or name" << std::endl
                  << "    Scale_factor is -1 for Auto-Scaling" << std::endl
                  << "    SSB_LSB is 1 for Upper Sideband, 0 for Lower Sideband" << std::endl;
        std::cout << std::endl; //blank line           
        // print the list of WaveOut devices
        const size_t numDevs = wave.getNumDevices();
        std::cout << "Found " << numDevs << " wave out devices." << std::endl;

        for (size_t k = 0; k < numDevs; ++k) {
            const std::string desc = wave.getDeviceDescription(k);
            std::cout << desc << std::endl;
        }
        return EXIT_FAILURE;
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
        std::cout << "Bad SSB Freq" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Using Receiver Band Index = " << nMem << std::endl;

    if ((sscanf(argv[2], "%d", &nWO) != 1)) {
        const std::string waveOutName(argv[2]);
        if (!wave.getDeviceByName(waveOutName, nWO)) {
            std::cout << "Bad WaveOut specified: " << waveOutName << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::cout << "Using wave output device: " << wave.getDeviceDescription(nWO) << std::endl;

    // Scale Factor
    if (argc >= 5) {
        if ((sscanf(argv[4], "%d", &SF) != 1) || (SF > 24)) {
            fprintf(stderr, "Bad Scale_factor (%s)\n", argv[4]);
            return EXIT_FAILURE;
        }
    }
    else {
        // Use AutoAF
        SF = -1;
    }
    if (-1 == SF) {
        std::cout << "Using automatic scale factor" << std::endl;
    }
    else {
        std::cout << "Using scale factor: " << SF << std::endl;
    }

    // USB/LSB.  USB = 1, LSB = 0
    int USB = 1;
    if ((sscanf(argv[3], "%d", &USB) != 1)) {
        fprintf(stderr, "Bad USB_LSB - must be 0 or 1 (%s)\n", argv[3]);
        return EXIT_FAILURE;
    }
    if (USB) {
        std::cout << "Demodulating Upper Sideband" << std::endl;
    }
    else {
        std::cout << "Demodulating Lower Sideband" << std::endl;
    }


    // Create AutoAF. Only used if SF == -1
    AutoScaleAF<float> af(2, wave.getClipValue());

    float scaleFactor = 0;
    float lastScaleFactor = 0;
    if (SF >= 0) {
        // User defined scale factor
        scaleFactor = pow(2, SF);
    }

    // create name of shared memory
    const std::string name = gPreSM + std::to_string(nMem) + gPostSM;

    // try to open shared memory
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

    // Create the demodulator

    const size_t SSB_BW = 3000;
    std::cout << "SSB Bandwidth: " << SSB_BW << " Hz" << std::endl;
    // F is always Fc-LO
    const int64_t LO = static_cast<int64_t>(SHDR->L0);
    const int64_t F = ssbFreq - SHDR->L0;
    SSBD<float> ssbd(SR, SSB_BW, static_cast<float>(F), static_cast<bool>(USB));
    const size_t SSB_SR = ssbd.GetOutRate();
    const size_t decRatio = SR / SSB_SR;
    const size_t ssbd_in_size = ssbd.GetInSize();

    // try to open waveout device
    const size_t Wave_SR = 48000;
    std::cout << "WaveOut Sample Rate=" << Wave_SR << std::endl;
    std::cout << "WaveOut Output Size=" << ssbd.GetOutSize() << std::endl;

    const float WAVE_BUFFER_TIME = 0.050;
    const size_t WAVE_BUFFER_LEN = static_cast<float>(Wave_SR) * WAVE_BUFFER_TIME; // WAVE_BUFFER_LEN is in samples
    std::cout << "Wave buffer length is " << WAVE_BUFFER_TIME << " seconds (" << WAVE_BUFFER_LEN << " samples)" << std::endl;
    const bool waveInitialized = wave.initialize(nWO, Wave_SR, WAVE_BUFFER_LEN);
    if (!waveInitialized) {
        std::cout << "Failed to open and start wave device." << std::endl;
        SM.Close();
        return EXIT_FAILURE;
    }

    Upsampler<float> upsamp(3); //2^3=8 ratio
    const size_t upsamp_ratio = upsamp.GetRatio();
    std::cout << "Upsampling wave audio from " << SSB_SR << " Hz to " << Wave_SR << " Hz" << std::endl;
    if (upsamp_ratio != Wave_SR/ SSB_SR) {
        std::cout << "Unsupported wave output sample rate" << std::endl;
        SM.Close();
        return EXIT_FAILURE;
    }

    std::cout << std::endl << "Press Q to terminate" << std::endl;

    const size_t iq_len = BIS;
    float* iq_data_f = reinterpret_cast<float*>(malloc(iq_len * 2 * sizeof(float)));

    //
    //  Main Loop
    //

    std::vector<float> af6khz(iq_len / decRatio, 0);
    std::vector<float> af48khz(iq_len * upsamp_ratio / decRatio, 0.0);

    std::vector<float> samples;
    samples.reserve(WAVE_BUFFER_LEN); // vector might grow beyond this

    bool terminate = false;

    while (!terminate) {
        // wait for new data. Blocks until data received
        SM.WaitForNewData();

        // read block of data from receiver

        const bool readSuccess = SM.Read((PBYTE)iq_data_f, iq_len * 2 * sizeof(float));
        if (!readSuccess) {
            std::cout << "Did not read any I/Q data" << std::endl;
            continue;
        }

        // Create complex from IQ
        std::complex<float> *xc = (std::complex<float>*)(iq_data_f);

        // Demodulate IQ for SSB
        for (size_t n = 0; n < iq_len; n += ssbd_in_size) {
            ssbd.Iterate(xc + n, af6khz.data() + n / decRatio);
        }

        // Upsample 6kHz audio to 48khz
        for (size_t n = 0; n < af6khz.size(); ++n) {
            upsamp.Iterate(af6khz.data() + n, af48khz.data() + n * upsamp_ratio);
        }

        // buffer the 48khz audio
        const size_t numWaveSamp = af48khz.size();
        for (size_t n = 0; n < numWaveSamp; ++n) {
            samples.push_back(af48khz[n]);
        }

        size_t numBufferedSamples = samples.size();
        // Write audio samples to waveout when we have enough
        if (numBufferedSamples >= WAVE_BUFFER_LEN) {
            size_t numSamplesToWrite = numBufferedSamples;
            // make sure we don't write more than WAVE_BUFFER_LEN samples
            if (numBufferedSamples > WAVE_BUFFER_LEN) {
                numSamplesToWrite = WAVE_BUFFER_LEN;
            }

            // If auto AF is enabled
            if (-1 == SF) {
                scaleFactor = af.getScaleFactor(samples);
                if (scaleFactor != lastScaleFactor) {
                    std::cout << "Scale factor changed to " << scaleFactor << std::endl;
                    lastScaleFactor = scaleFactor;
                }
            }

            const bool writeSuccess = wave.write(samples.data(), numSamplesToWrite, scaleFactor);

            // Remove written samples
            if (numSamplesToWrite == numBufferedSamples) {
                // Entire vector was written, so clear all
                samples.clear();
            }
            else {
                // Partial vector write. Remove written samples only
                samples.erase(samples.begin(), samples.begin() + numSamplesToWrite);
            }
        }

        // Was Exit requested?
        if (_kbhit()) {
            const char ch = _getch();
            if (ch == 'Q' || ch == 'q') {
                std::cout << "Q pressed, so terminating" << std::endl;
                terminate = true;
            }
        }
    }

    return EXIT_SUCCESS;

    wave.stop();
    SM.Close();
    free(iq_data_f);
}
