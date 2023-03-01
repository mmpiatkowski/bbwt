#ifndef _BBWT_HPP_
#define _BBWT_HPP_

/*
 * Bijective Burrows-Wheeler Transform implementation.
 *
 * (c) 2023 Marcin Piątkowski, marcin.piatkowski(at)mat.umk.pl
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "BitVector.hpp"
#include "lyndon.hpp"
#include "bbwt_internal.hpp"


/** Computes the circular suffix array of inStr.
 * @param inStr input data buffer
 * @param csa buffer where computed circular suffix array is stored
 * @param len the size of the input data
 * @param alphSize size of the alphabet
 * @return 0 after successful computation, non-zero in case of any error
 */
template<typename Tdata, typename Tnum>
int circularSuffixArray(const Tdata *inStr, Tnum *csa, Tnum len, const Tnum alphSize = 256) {

    //------------------------------------------------------------------------------------------------------------------
    // Incorrect and trivial input data
    //------------------------------------------------------------------------------------------------------------------

    if (inStr == nullptr || csa == nullptr) {
        return -1;
    }

    if (len == 0)
        return 0;

    if (len == 1) {
        csa[0] = 0;

        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Compute Lyndon factorisation of the input data
    //------------------------------------------------------------------------------------------------------------------

    BitVector lFac(len + 1);
    lyndonFactors(inStr, len, &lFac);

    //------------------------------------------------------------------------------------------------------------------
    // Compute Circular Suffix Array using modified SAIS algorithm
    //------------------------------------------------------------------------------------------------------------------

    return circularSuffixArray(inStr, csa, len, lFac, alphSize);
}



/** Computes Bijective Burows-Wheeler Transform of inStr.
 * @param inStr input data buffer
 * @param outStr buffer where the computed BBWT is stored (may be the same as inStr)
 * @param csa memory buffer where circular suffix array will be stored
 * @param len the size of the input data
 * @param alphSize size of the alphabet
 * @return 0 after successful computation, non-zero in case of any error
 */
template<typename Tdata, typename Tnum>
int bbwt(const Tdata *inStr, Tdata *outStr, Tnum *csa, Tnum len, const Tnum alphSize = 256) {

    //------------------------------------------------------------------------------------------------------------------
    // Incorrect and trivial input data
    //------------------------------------------------------------------------------------------------------------------

    if (inStr == nullptr || outStr == nullptr) {
        return -1;
    }

    if (len == 0)
        return 0;

    if (len == 1) {
        outStr[0] = inStr[0];

        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Compute Lyndon factorisation of the input data
    //------------------------------------------------------------------------------------------------------------------

    BitVector lFac(len + 1);    // All Lyndon factors
    BitVector lFirst(len + 1);  // Only the first occurrence of each Lyndon factor

    lyndonFactors(inStr, len, &lFac, &lFirst);

    //------------------------------------------------------------------------------------------------------------------
    // Compute circular suffix array for the input data
    //------------------------------------------------------------------------------------------------------------------

    circularSuffixArray(inStr, csa, len, lFac, alphSize);

    //------------------------------------------------------------------------------------------------------------------
    // Retrieve Bijective Burrows-Wheeler Transform from circular suffix array.
    // If the input and output buffer have a non-empty overlap use suffix array buffer as a temporary storage.
    //------------------------------------------------------------------------------------------------------------------

    if (inStr > outStr + len || outStr > inStr + len) {
        for (Tnum outPos = 0; outPos < len; ++outPos) {
            Tnum inPos = csa[outPos];

            // Wrap around the Lyndon factor if needded
            if (lFac.get(inPos)) {
                inPos = lFac.next(inPos) - 1;
            } 
            else {
                --inPos;
            }

            outStr[outPos] = inStr[inPos];
        }
    }
    else {
        for (Tnum outPos = 0; outPos < len; ++outPos) {
            Tnum inPos = csa[outPos];

            // Wrap around the Lyndon factor if needed
            if (lFac.get(inPos)) {
                inPos = lFac.next(inPos) - 1;
            } 
            else {
                --inPos;
            }

            csa[outPos] = inStr[inPos];
        }

        for (Tnum pos = 0; pos < len; ++pos) {
            outStr[pos] = (unsigned char) csa[pos];
        }
    }

    return 0;
}


/** Computes the inverse of Bijective Burrows-Wheeler Transform of inStr.
 * @param inStr input data
 * @param outStr buffer where the computed inverse of BBWT is stored
 * @param len the size of the input data
 * @param alphSize size of the alphabet
 * @return 0 after successful computation, non-zero in case of any error
 */
template<typename Tdata, typename Tnum>
int unbbwt(const Tdata *inStr, Tdata *outStr, Tnum len, const Tnum alphSize = 256) {
    const Tnum MaxVal = std::numeric_limits<Tnum>::max();

    //------------------------------------------------------------------------------------------------------------------
    // Incorrect and trivial input data
    //------------------------------------------------------------------------------------------------------------------

    if (inStr == nullptr || outStr == nullptr) {
        return -1;
    }

    if (len == 0)
        return 0;

    if (len == 1) {
        outStr[0] = inStr[0];

        return 0;
    }

    Tnum charsCount[alphSize] = {0};
    Tnum charsBefore[alphSize] = {0};
    Tnum charsSeen[alphSize] = {0};
    Tnum *stdPerm = new Tnum[len];

    for (Tnum i=0; i<len; ++i)
        ++charsCount[inStr[i]];

    for (Tnum i=1; i<alphSize; ++i)
        charsBefore[i] = charsBefore[i-1] + charsCount[i-1];

    for (Tnum i=0; i<len; ++i) {
        stdPerm[i] = charsBefore[inStr[i]] + charsSeen[inStr[i]];
        ++charsSeen[inStr[i]];
    }

    Tnum outPos = len - 1;

    for (Tnum j=0; j<len; ++j) {
        if (stdPerm[j] == MaxVal)
            continue;

        Tnum inPos = j;

        while (stdPerm[inPos] != MaxVal) {
            outStr[outPos] = inStr[inPos];
            --outPos;
            Tnum t = inPos;
            inPos = stdPerm[inPos];
            stdPerm[t] = MaxVal;
        }
    }

    delete[] stdPerm;

    return 0;
}


#endif //_BBWT_HPP_