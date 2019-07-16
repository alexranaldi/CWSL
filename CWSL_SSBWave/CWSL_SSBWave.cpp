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

std::atomic_bool terminateFlag;
std::atomic_bool holdScaleFactor;

// Used internally as a single producer single consumer queue (ring buffer)
template <typename T>
struct ring_buffer_t {
    T* recs;
    // Use atomics for indices to prevent instruction reordering
    std::atomic_int read_index;
    std::atomic_int write_index;
    // Number of items the ring buffer can hold.
    size_t size;
};

template <typename T>
bool initialize(ring_buffer_t<T>& buf, const size_t size) {
    buf.read_index = 0;
    buf.write_index = 0;
    buf.size = size;
    buf.recs = reinterpret_cast<T*>(malloc(sizeof(T) * buf.size));
    return (buf.recs != nullptr);
}

template <typename T>
bool wait_for_empty_slot(ring_buffer_t<T>& buf) {
    while ((buf.read_index == buf.write_index + 1) || (buf.read_index == 0 && static_cast<int64_t>(buf.write_index) == static_cast<int64_t>(buf.size) - 1)) {
        std::this_thread::sleep_for(std::chrono::microseconds(100));
        if (terminateFlag) {
            return false;
        }
    }
    return true;
}

template <typename T>
void inc_write_index(ring_buffer_t<T>& buf) {
    // cast to signed so we don't break subtraction
    if (static_cast<int64_t>(buf.write_index) == static_cast<int64_t>(buf.size) - 1) {
        buf.write_index = 0;
    }
    else {
        buf.write_index++;
    }
}

template <typename T>
T pop(ring_buffer_t<T>& buf) {
    wait_for_data(buf);
    T curr = buf.recs[buf.read_index];
    if (buf.read_index == buf.size - 1) {
        buf.read_index = 0;
    }
    else {
        buf.read_index++;
    }
    return curr;
}

template <typename T>
bool wait_for_data(ring_buffer_t<T>& buf) {
    while (buf.read_index == buf.write_index) {
        std::this_thread::sleep_for(std::chrono::microseconds(100));
        if (terminateFlag) {
            return false;
        }
    }
    return true;
}

ring_buffer_t<std::complex<float>*> iq_buffer;
ring_buffer_t<float*> af_buffer;

int SF = 16;

int SMNumber = -1;

void readIQ(CSharedMemory &SM, const size_t iq_len) {

    while (!terminateFlag) {
    
        // wait for new data from receiver. Blocks until data received
        SM.WaitForNewData();
        
        // wait for a slot to be ready in the buffer. Blocks until slot available.
        if (!wait_for_empty_slot(iq_buffer)) {
            std::cout << "No slots available in IQ buffer!" << std::endl;
            continue;
        }

        std::complex<float>* iq_data_f = iq_buffer.recs[iq_buffer.write_index];

        // read block of data from receiver
        const bool readSuccess = SM.Read((PBYTE)iq_data_f, iq_len * sizeof(std::complex<float>));
        if (readSuccess) {
            inc_write_index(iq_buffer);
        }
        else {
            std::cout << "Did not read any I/Q data from shared memory" << std::endl;
        }
    }
}

void demodulate(SSBD<float>& ssbd, Upsampler<float>& upsamp, AutoScaleAF<float>& af, const size_t iq_len, const size_t decRatio) {
    // Demodulate IQ for SSB

    const size_t ssbd_in_size = ssbd.GetInSize();
    const size_t upsamp_ratio = upsamp.GetRatio();

    std::vector<float> af6khz(iq_len / decRatio, 0.0);

    const size_t n48khz = iq_len * upsamp_ratio / decRatio;

    float scaleFactor = static_cast<float>(pow(SF, 2));
    float lastScaleFactor = 0;

    while (!terminateFlag) {
        // wait for available memory to write audio data to
        if (!wait_for_empty_slot(af_buffer)) {
            std::cout << "No slots available for audio data!" << std::endl;
            continue;
        }
        // get IQ data
        std::complex<float> *xc = pop(iq_buffer);
        // The audio data we're going to write to
        float* af48khz = af_buffer.recs[af_buffer.write_index];

        for (size_t n = 0; n < iq_len; n += ssbd_in_size) {
            ssbd.Iterate(xc + n, af6khz.data() + n / decRatio);
        }

        // Upsample 6kHz audio to 48khz
        for (size_t n = 0; n < af6khz.size(); ++n) {
            upsamp.Iterate(af6khz.data() + n, af48khz + n * upsamp_ratio);
        }

        // If auto AF is enabled
        if (-1 == SF && !holdScaleFactor) {
            std::vector<float> afSamples(af48khz, af48khz + n48khz);
            scaleFactor = af.getScaleFactor(afSamples);
            if (scaleFactor != lastScaleFactor) {
                std::cout << "Scale factor set to: " << scaleFactor << std::endl;
                lastScaleFactor = scaleFactor;
            }
        }
        for (size_t n = 0; n < n48khz; ++n) {
            af48khz[n] *= scaleFactor;
        }
        // AF write complete
        inc_write_index(af_buffer);

    }
}

void writeAudio(WinWave& wave, const size_t len) {

    while (!terminateFlag){
        // get af data
        float* af = pop(af_buffer);
        // write data to wave device
        const bool writeSuccess = wave.write(af, len);
        if (!writeSuccess) {
            std::cout << "Failed to write to audio device" << std::endl;
        }
    }

}

std::string createSharedMemName(const int bandIndex) {
    // create name of shared memory
    std::string Name = gPreSM + std::to_string(bandIndex) + gPostSM;
    if (SMNumber != -1) {
        Name += std::to_string(SMNumber);
    }
    return Name;
}

//////////////////////////////////////////////////////////////////////////////
// Find the right band
int findBand(const int64_t f) {
    CSharedMemory SM;
    SM_HDR h;

    // try to find right band - for all possible bands ...
    for (int bandIndex = 0; bandIndex < MAX_CWSL; ++bandIndex) {

        const std::string Name = createSharedMemName(bandIndex);

        std::cout << "Opening shared memory with name: " << Name << std::endl;

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
    int nMem = 0;
    size_t nWO = 0;

    terminateFlag = false;

    WinWave wave = WinWave();

    // check number of input parameters
    if (argc < 4) {
        // print usage
        std::cout << "Not enough input arguments!" << std::endl;
        std::cout << "Usage: CWSL_SSBWave FreqHz WaveOutNr SSB_LSB Scale_factor Shared_Mem" << std::endl;
        std::cout << "    FreqHz is the frequency in Hz" << std::endl
                  << "    WaveOutNr is the Wave Out device number or name" << std::endl
                  << "    SSB_LSB is 1 for Upper Sideband, 0 for Lower Sideband" << std::endl
                  << "    Scale_factor is -1 for Auto-Scaling" << std::endl
                  << "    Shared_Mem is an optional single numeric digit specifying the shared memory interface" << std::endl;
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

    // Shared Mem
    if (argc >= 6) {
        std::cout << "A shared memory interface was specified." << std::endl;
        if (sscanf(argv[5], "%d", &SMNumber) != 1) {
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
        std::cout << "Using automatic scale factor (AutoAF)" << std::endl;
    }
    else {
        std::cout << "Using user-specified scale factor: " << SF << std::endl;
    }

    // USB/LSB.  USB = 1, LSB = 0
    int USB = 1;
    if ((sscanf(argv[3], "%d", &USB) != 1)) {
        fprintf(stderr, "Bad USB_LSB - must be 0 or 1 (%s)\n", argv[3]);
        return EXIT_FAILURE;
    }
    if (USB) {
        std::cout << "Demodulating Upper Sideband." << std::endl;
    }
    else {
        std::cout << "Demodulating Lower Sideband." << std::endl;
    }

    //
    // Setup shared memory and receiver interface
    //
   
    // create name of shared memory
    const std::string name = createSharedMemName(nMem);
    // try to open shared memory
    if (!SM.Open(name.c_str())) {
        fprintf(stderr, "Can't open shared memory for %d receiver\n", nMem);
        std::cout << "Failed to open shared memory, so terminating" << std::endl;
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
    const bool waveInitialized = wave.initialize(nWO, Wave_SR, af_len);
    if (!waveInitialized) {
        std::cout << "Failed to open and start wave device." << std::endl;
        SM.Close();
        return EXIT_FAILURE;
    }
    std::cout << "Opened Audio Device" << std::endl;
    std::cout << "Audio Device Sample Rate=" << Wave_SR << std::endl;
    const auto waveBPS = wave.getBitsPerSample();
    std::cout << "Audio Device Bits Per Sample=" << waveBPS << std::endl;
    wave.enablePrintClipped();


    //
    // Create AutoAF. Only used if SF == -1
    //

    const float headroomdB = 18;
    const float windowdB = 22;
    const size_t numSampBeforeAfInc = Wave_SR * 300; // 300s of samples
    const float clipVal = wave.getClipValue();
    std::cout << "Wave device absolute maximum signal value (clip level): " << clipVal << std::endl;
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

    initialize(iq_buffer, 4096);
    for (size_t k = 0; k < iq_buffer.size; ++k) {
        iq_buffer.recs[k] = reinterpret_cast<std::complex<float>*>( malloc( sizeof(std::complex<float>) * iq_len ) );
    }
    initialize(af_buffer, 4096);
    for (size_t k = 0; k < af_buffer.size; ++k) {
        af_buffer.recs[k] = reinterpret_cast<float*>( malloc( sizeof(float) * af_len) );
    }

    //
    //  Start Threads
    //

    std::cout << "Creating receiver thread..." << std::endl;
    std::thread iqThread = std::thread(readIQ, std::ref(SM), BIS);
    std::cout << "Creating SSB Demodulator thread..." << std::endl;
    std::thread demodThread = std::thread(&demodulate, std::ref(ssbd), std::ref(upsamp), std::ref(af), BIS, decRatio);
    std::cout << "Creating Audio Device thread..." << std::endl;
    std::thread audioThread = std::thread(&writeAudio, std::ref(wave), af_len);
    
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
    
    if (audioThread.joinable()) {
        audioThread.join();
    }
    if (demodThread.joinable()) {
        demodThread.join();
    }
    if (iqThread.joinable()) {
        iqThread.join();
    }
    


    //
    //  Clean up
    //

    wave.stop();
    SM.Close();

    std::cout << "Deallocating memory..." << std::endl;

    // free all allocated memory
    for (size_t k = 0; k < af_buffer.size; ++k) {
        free(af_buffer.recs[k]);
    }
    for (size_t k = 0; k < iq_buffer.size; ++k) {
        free(iq_buffer.recs[k]);
    }
    free(af_buffer.recs);
    free(iq_buffer.recs);

    return EXIT_SUCCESS;

}    
