// Alex Ranaldi  W2AXR   alexranaldi@gmail.com

// LICENSE: GNU General Public License v3.0
// THE SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.

#include <string>
#include <iostream>
#include <limits>
#include <cmath>

#include <Windows.h>

#include "WinWave.hpp"

///////////////////////////////////////////////////////////////////////////////
// Construction
WinWave::WinWave() : 
    mInitialized(false),
    mFs(0),
    mBitsPerSample(0),
    mClipValue(0),
    mPrintClipped(false) {
    mBuffers.resize(NUM_BUFFERS);
    mHeaders.resize(NUM_BUFFERS);
}

WinWave::~WinWave() {
    stop();
}

bool WinWave::initialize(const size_t devNum, const uint64_t Fs, const size_t bufferLen) {
    // if already initialized
    if (mInitialized) {
        return false;
    }

    bool success = false;
    size_t bitsPerSample = 32;
    while ( (!success) && (bitsPerSample >= 8) ){
        success = openWaveDevice(devNum, Fs, bitsPerSample);
        if (!success) {
            bitsPerSample -= 8;
        }
    }
        
    // Was the device initialized?
    if (success) {
        allocBuffers(bufferLen);
        mInitialized = true;
        mClipValue = static_cast<float>(std::pow(2, bitsPerSample - 1) - 1);
        return true;
    }
    
    // Device not opened
    return false;
}

size_t WinWave::getBitsPerSample() const{
    return mBitsPerSample;
}

void WinWave::allocBuffers(const size_t bufferLen) {
    mBufferLen = bufferLen;
    for (size_t k = 0; k < NUM_BUFFERS; ++k) {
        mBuffers[k] = (int32_t*)malloc(sizeof(int32_t) * mBufferLen);
    }
    mBufferIndex = 0;
}

bool WinWave::openWaveDevice(const size_t devNum, const uint64_t Fs, const size_t bitsPerSample) {
    bool initialized = false;
   
    WAVEFORMATEX format;

    format.wFormatTag = WAVE_FORMAT_PCM;
    format.nChannels = NUM_CHANNELS;
    format.nSamplesPerSec = static_cast<DWORD>(Fs);
    format.wBitsPerSample = static_cast<WORD>(bitsPerSample);
    // nBlockAlign = num chans * bits/sample / 8, per Microsoft doc
    format.nBlockAlign = format.nChannels * format.wBitsPerSample / 8;
    format.nAvgBytesPerSec = format.nSamplesPerSec * format.nBlockAlign;
    format.cbSize = 0; // ignored

    // Open the Wave Out audio device
    const MMRESULT res = waveOutOpen(&mWaveOut, devNum, &format, NULL, NULL, CALLBACK_NULL);
    initialized = res == MMSYSERR_NOERROR;
    if (initialized) {
        mFs = Fs;
        mDevNum = devNum;
        mBitsPerSample = bitsPerSample;
    }
    return initialized;
}

float WinWave::getClipValue() const{
    return mClipValue;
}

void WinWave::enablePrintClipped(){
    mPrintClipped = true;
}

bool WinWave::write(float* samples, const size_t numSamples) {

    if (numSamples > mBufferLen) {
        return false;
    }

    const size_t nChannels = 1;

    WAVEHDR& header = mHeaders[mBufferIndex];
    int32_t* buffer = mBuffers[mBufferIndex];

    for (size_t k = 0; k < numSamples; ++k) {
        const float val = samples[k];
        if ( mPrintClipped && ((val >= mClipValue) || (val <= -mClipValue)) ) {
            std::cout << "Clip" << std::endl;
        }
        // float -> int32_t
        buffer[k] = static_cast<int32_t>(val);
    }

    header.lpData = (LPSTR)buffer;
    header.dwBufferLength = numSamples * sizeof(int32_t);
    header.dwFlags = 0;

    MMRESULT mmr = waveOutPrepareHeader(mWaveOut, &header, sizeof(WAVEHDR));
    bool ok = mmr == MMSYSERR_NOERROR;
    if (!ok) {
        return false;
    }
    mmr = waveOutWrite(mWaveOut, &header, sizeof(WAVEHDR));
    ok = mmr == MMSYSERR_NOERROR;
    if (ok) {
        // Go to next buffer
        mBufferIndex++;
        if (mBufferIndex == NUM_BUFFERS - 1) {
            mBufferIndex = 0;
        }
    }

    return ok;
}

std::string WinWave::getDeviceDescription(const size_t devNum) {
    WAVEOUTCAPS caps = {};
    const MMRESULT mmr = waveOutGetDevCaps(devNum, &caps, sizeof(caps));
    if (MMSYSERR_NOERROR != mmr) {
        return "";
    }
    return std::to_string (devNum) + " " + caps.szPname;
}

std::string WinWave::getDeviceName(const size_t devNum) {
    WAVEOUTCAPS caps = {};
    const MMRESULT mmr = waveOutGetDevCaps(devNum, &caps, sizeof(caps));
    if (MMSYSERR_NOERROR != mmr) {
        return "";
    }
    return std::string(caps.szPname);
}

bool WinWave::getDeviceByName(const std::string devName, size_t &devNum) {
    const size_t numDevices = getNumDevices();
    for (size_t k = 0; k < numDevices; ++k) {
        const std::string testName = getDeviceName(k);
        const int same = testName.compare(0, devName.length(), devName);
        // 0 means strings equal
        if (0 == same) {
            devNum = k;
            return true;
        }
    }
    return false;
}

size_t WinWave::getNumDevices() {
    return waveOutGetNumDevs();
}

void WinWave::stop() {
    if (mInitialized){
        mInitialized = false;
        // Free memory in each buffer
        for (size_t k = 0; k < NUM_BUFFERS; ++k) {
            free(mBuffers[k]);
        }
    }
}
