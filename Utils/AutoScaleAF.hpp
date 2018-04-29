// Alex Ranaldi  W2AXR   alexranaldi@gmail.com

// LICENSE: GNU General Public License v3.0
// THE SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.

#include <vector>
#include <algorithm>

// Auto computes a factor that can be used to scale audio data.
//  This is essentially a volume knob that turns itself.
template <class T>
class AutoScaleAF{
    
    public:
            
        // headroomdB specifies the ideal headroom in dB to maintain
        // maxVal specifies the max. e.g., a value that represents a clipped value
        AutoScaleAF(const T headroomdB, const T headroomWindowdB, const T maxVal, const size_t adjSampInterval) :
            mMaxVal(maxVal),
            mMaxdB(0),
            mScaleFactor(0),
            mHiThreshdB(0),
            mLoThreshdB(0),
            mAmountDecrease(0),
            mAmountIncrease(0),
            mNumSampBelowMinThresh(0),
            mAdjSampInterval(adjSampInterval),
            mNumSampSinceDecrease(0) {
                
            // maximum value before signal clips
            mMaxdB = v_to_db(maxVal);
            mHiThreshdB = mMaxdB - headroomdB;
            mLoThreshdB = mHiThreshdB - headroomWindowdB;
            mHiThreshV = db_to_v(mHiThreshdB);
            mLoThreshV = db_to_v(mLoThreshdB);
            mAmountDecrease = 1/db_to_v(3);
            mAmountIncrease = db_to_v(1);
        }
        
        virtual ~AutoScaleAF(){
        }
    
        T getScaleFactor(const std::vector<T> v){
            computeScaleFactor(v);
            return mScaleFactor;
        }

        void getThresholds(T& hi, T& lo){
            hi = mMaxdB - mHiThreshdB;
            lo = mMaxdB - mLoThreshdB;
        }
        
    private:
    
        // accepts samples and recomputes the scale factor
        void computeScaleFactor(const std::vector<T> v){
            const size_t numSamp = v.size();
            // If no scale factor has been determined yet
            if (0 == mScaleFactor){
                // initial scale factor
                const auto maxIt = max_element(v.begin(), v.end());
                const T maxVal = *maxIt;
                mScaleFactor = mLoThreshV / maxVal;
            }
            else {
                const size_t numGTMax = num_greater_than_abs(v, mMaxVal / mScaleFactor);
                const size_t numGTHi = num_greater_than_abs(v, mHiThreshV / mScaleFactor);
                const T percentGTHi = static_cast<T>(numGTHi) / static_cast<T>(numSamp);

                // If more than 2.5% clip, or more than 10% of the samples are above the ideal max
                if (numGTMax >= 0.025 || percentGTHi >= 0.10){
                    // Decrease sale factor
                    mScaleFactor *= mAmountDecrease;
                    mNumSampSinceDecrease = 0;
                    mNumSampBelowMinThresh = 0;
                }
                else {
                    mNumSampSinceDecrease += numSamp;
                    if (mNumSampSinceDecrease >= mAdjSampInterval){
                        const size_t numGTLo = num_greater_than_abs(v, mLoThreshV / mScaleFactor);
                        const T percentGTLo = static_cast<T>(numGTLo) / static_cast<T>(numSamp);
                        const T percentLTLo = 1 - percentGTLo;
                        // If more than 50% of samples are below the ideal min
                        // (aka if less than 50% are above the ideal min) 
                        if (percentLTLo > 0.50) {
                            mNumSampBelowMinThresh += numSamp;
                        }
                        else {
                            mNumSampBelowMinThresh = 0;
                        }
                        if (mNumSampBelowMinThresh > mAdjSampInterval / 2) {
                            mScaleFactor *= mAmountIncrease;
                            mNumSampBelowMinThresh = 0;
                        }
                    }
                }
            }
        }
    
        T v_to_db(const T v_in) const{
            return 20*std::log10(v_in);
        }
        
        T db_to_v(const T db_in) const{
            return static_cast<T>(std::pow(10, db_in / 20));
        }

        size_t num_greater_than_abs(const std::vector<T> v, const T test) const{
            size_t count = 0;
            for (size_t k = 0; k < v.size(); ++k){
                if ( (v[k] > test) || (v[k] < -test) ) {
                    count++;
                }
            }
            return count;
        }

        T mMaxVal;
        T mMaxdB;
        T mScaleFactor;
        T mHiThreshV;
        T mLoThreshV;
        T mHiThreshdB;
        T mLoThreshdB;
        T mAmountDecrease;
        T mAmountIncrease;
        size_t mNumSampSinceDecrease;
        size_t mNumSampBelowMinThresh;
        size_t mAdjSampInterval;

 };