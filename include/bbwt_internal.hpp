#ifndef _BBWT_INTERNAL_HPP
#define _BBWT_INTERNAL_HPP

/**
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

#include <new>
#include <limits>
#include <vector>

#include "BitVector.hpp"
#include "lyndon.hpp"


//-------------------------------------------------------------------------------------------------
// Buckets
//-------------------------------------------------------------------------------------------------

/**
 * Computes computes begins and ends of all buckets related to inData.
 *
 * @param inData The buffer containing data for which the buckets structure should be computed
 * @param len The size of the input data.
 * @param buckets The output buckets structure. The bucket related to the i-th character
 *                (in the lexicographical order) is set as buckets[i]..buckets[i+1].
 *                The array should have no less than len+1 elements.
 * @param alphSize The size of the alphabet considered
 */
template<typename Tdata, typename Tnum>
void computeBucketsStructure(const Tdata *inData, const Tnum len, Tnum *buckets, const Tnum alphSize = 256) {
    std::vector<Tnum> charCount(alphSize, 0);
    Tnum total = 0;

    for (Tnum i=0; i<len; ++i) {
        ++charCount[inData[i]];
    }

    buckets[0] = 0;

    for (Tnum i=0; i<alphSize; ++i) {
        total += charCount[i];
        buckets[i+1] = total;
    }
}


//-------------------------------------------------------------------------------------------------
// Circular Suffix Array Construction
//-------------------------------------------------------------------------------------------------

// Types of suffixes.
const int LType = 0;
const int SType = 1;


/** Returns "true" if "pos" is a starTdatag position of LMS inf-substring and "false" otherwise. */
template<typename Tnum>
bool isLMSPos(const Tnum pos, const BitVector<Tnum> & lFac, const BitVector<Tnum> & sType) {
    return lFac.get(pos) || (sType.get(pos)==SType && sType.get(pos - 1) == LType);
}

/** Returns "true" if the Lyndon factor starTdatag at pos "pos" has length 1 and "false" otherwise. */
template<typename Tnum>
bool isSingleton(const Tnum pos, const BitVector<Tnum> & lFac) {
    return lFac.get(pos) && lFac.next(pos) == (pos + 1);
}


/*
 * Place all suffixes of type L at the beginning of corresponding bucket.
 */
template<typename Tdata, typename Tnum>
int preSortSuffixexL(const Tdata *inStr, Tnum *sa, Tnum len, const BitVector<Tnum> &lFac, const BitVector<Tnum> &suffType, const BitVector<Tnum> &spcFac, Tnum *buckets, const Tnum alphSize=256) {
    Tnum p = spcFac.prev(len);

    for (Tnum i=0; i<len; ++i) {
        while ( (p >= 0) && buckets[inStr[p]] == i) {
            Tnum j = lFac.next(p) - 1;
            sa[buckets[inStr[j]]] = j;
            ++buckets[inStr[j]];
            p = spcFac.prev(p);
        }

        Tnum j = sa[i];

        if (j < 0) {
            continue;
        }

        // Wrap around Lyndon factor if needed
        if (lFac.get(j) == 0) {
            --j;
        }
        else {
            j = lFac.next(j) - 1;
        }

        if (suffType.get(j) == LType) {
            sa[buckets[inStr[j]]] = j;
            ++buckets[inStr[j]];
        }
    }

    return 0;
}


/*
 * Place all suffixes of type S at the end of corresponding bucket.
 */
template<typename Tdata, typename Tnum>
int preSortSuffixesS(const Tdata *inStr, Tnum *sa, Tnum len, const BitVector<Tnum> &lFac, const BitVector<Tnum> &suffType, Tnum *buckets, const Tnum alphSize=256) {
    for (Tnum i=len-1; i>=0; --i) {
        Tnum j = sa[i];

        if (j < 0) {
            continue;
        }

        if (lFac.get(j) == 0) {
            --j;

            if (suffType.get(j) == SType) {
                --buckets[inStr[j] + 1];
                sa[buckets[inStr[j] + 1]] = j;
            }
        }
    }

    return 0;
}


template<typename Tdata, typename Tnum>
int circularSuffixArray(const Tdata *inStr, Tnum *sa, Tnum len, const BitVector<Tnum> &lbFac, const Tnum alphSize = 256) {

    //------------------------------------------------------------------------------------------------------------------
    // Mark each position (and corresponding suffix) in inStr as type S or L respectively.
    // All suffixes are initially assumed to be of type L (0), therefore we need to mark type S (1) suffixes only.
    //------------------------------------------------------------------------------------------------------------------

    BitVector<Tnum> suffType(len + 7);
    BitVector<Tnum> spcSuff(len + 1);

    for (Tnum fStart=0, fEnd; fStart<len; fStart=fEnd) {
        fEnd = lbFac.next(fStart);
        suffType.set(fStart, SType);

        // Mark all suffixes of type S
        for (Tnum j = fEnd - 2; j >= fStart; --j) {
            if ((inStr[j] < inStr[j + 1]) || (inStr[j] == inStr[j + 1] && suffType.get(j + 1) == SType))
                suffType.set(j, SType);
        }

        Tnum m = 0, c = 0, c0;
        Tnum c1 = inStr[fEnd - 1];

        for (Tnum i = fEnd - 2; i >= fStart; --i) {
            if ((c0 = inStr[i]) < (c1 + c)) {
                c = 1;
            } 
            else if (c != 0) {
                m += 1;
                c = 0;
            }

            c1 = c0;
        }

        if ((m == 0) && (c == 0)) {
            spcSuff.set(fStart, true);
        }
    }

    spcSuff.set(len, true);

    //------------------------------------------------------------------------------------------------------------------
    // Compute bucket sizes for the input data
    //------------------------------------------------------------------------------------------------------------------

    Tnum *buckets;
    Tnum *tmpBuckets;

    try {
        buckets = new Tnum[alphSize + 1];
        tmpBuckets = new Tnum[alphSize + 1];
    }
    catch (const std::bad_alloc &e) {
        return -1;
    }

    computeBucketsStructure(inStr, len, buckets, alphSize);

    // Initialise all suffixes as being not set in a proper order
    for (Tnum i=0; i<len; ++i) {
        sa[i] = -1;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Insert each LMS inf-suffix (except those of length 1) at the end of the corresponding bucket
    //------------------------------------------------------------------------------------------------------------------

    memcpy(tmpBuckets, buckets, (alphSize+1)*sizeof(Tnum));

    for (Tnum i=0; i<len; ++i) {
        if (isLMSPos(i, lbFac, suffType) && !spcSuff.get(i)) {
            sa[tmpBuckets[inStr[i]+1]-1] = i;
            --tmpBuckets[inStr[i]+1];
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    // Insert L inf-suffixes into the proper bucket (starting from the beginning of the bucket)
    //------------------------------------------------------------------------------------------------------------------
    memcpy(tmpBuckets, buckets, (alphSize+1)*sizeof(Tnum));
    preSortSuffixexL(inStr, sa, len, lbFac, suffType, spcSuff, tmpBuckets, alphSize);

    //------------------------------------------------------------------------------------------------------------------
    // Insert S inf-suffixes into the proper bucket (starting from the bucket end)
    //------------------------------------------------------------------------------------------------------------------
    memcpy(tmpBuckets, buckets, (alphSize+1)*sizeof(Tnum));
    preSortSuffixesS(inStr, sa, len, lbFac, suffType, tmpBuckets, alphSize);

    //------------------------------------------------------------------------------------------------------------------
    // Compact all LMS inf-suffixes into the first positions in the suffix array and clear its remaining part
    // To reduce the space complexity we use the end of the suffix array buffer to store labels of LMS inf-suffixes
    //------------------------------------------------------------------------------------------------------------------

    Tnum numLMSSuff = 0;

    for (Tnum i=0; i<len; ++i) {
        if (isLMSPos(sa[i], lbFac, suffType) && !spcSuff.get(sa[i])) {
            sa[numLMSSuff] = sa[i];
            ++numLMSSuff;
        }
    }

    for (Tnum i=numLMSSuff; i<len; ++i)
        sa[i] = 0;

    //------------------------------------------------------------------------------------------------------------------
    // Compute meta-labels for LMS inf-suffixes preserving lexicographic order
    // Encode LMS inf-suffixes with the proper meta-labels into the shorter string
    // (of the length at most half the length of the original string)
    //------------------------------------------------------------------------------------------------------------------

    Tnum *labelBuff;
    try {
        labelBuff = new Tnum[numLMSSuff + 1];
    }
    catch (const std::bad_alloc &e) {
        return -1;
    }

    for (Tnum fStart=0, fEnd; fStart<len; fStart=fEnd) {
        fEnd = lbFac.next(fStart);

        Tnum j = fEnd, c = 0, c0 = 0, c1 = inStr[fEnd - 1];

        for (Tnum i = fEnd - 2; i >= fStart; --i) {
            if ((c0 = inStr[i]) < (c1 + c)) {
                c = 1;
            }
            else if(c != 0) {
                sa[numLMSSuff + ((i + 1) >> 1)] = j - i - 1;
                j = i + 1;
                c = 0;
            }

            c1 = c0;
        }

        if ((j < fEnd) || (c != 0)) {
            sa[numLMSSuff + (fStart >> 1)] = j - fStart;
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    // Find the meta-labels of all LMS inf-suffixes preserving their lexicographical order
    //------------------------------------------------------------------------------------------------------------------

    Tnum numLabels = 0;
    Tnum q = len;
    Tnum qLen = 0;

    for (Tnum i=0; i<numLMSSuff; ++i) {
        Tnum pos = sa[i];
        Tnum sbwrdLen = sa[numLMSSuff + (pos >> 1)];
        bool distinct = true;

        if (sbwrdLen == qLen) {
            Tnum j = 0;

            for (j = 0; (j < sbwrdLen) && (inStr[pos + j] == inStr[q + j]); ++j) { }

            if(j == sbwrdLen) {
                distinct = false;
            }
        }

        if (distinct) {
            ++numLabels;
            q = pos;
            qLen = sbwrdLen;
        }

        sa[numLMSSuff + (pos >> 1)] = numLabels;
    }

    //------------------------------------------------------------------------------------------------------------------
    // If the labels of LMS inf-suffixes are not unique we need a recursive call to sort them properly
    //------------------------------------------------------------------------------------------------------------------

    if (numLabels < numLMSSuff) {
        // Derive the Lyndon factorisation of the reduced string from the Lyndon factorisation of the original string
        BitVector<Tnum> redFactors(numLMSSuff + 1);

        for (Tnum inPos=0, outPos=0; inPos < len; ++inPos) {
            if (isLMSPos(inPos, lbFac, suffType) && !spcSuff.get(inPos)) {
                redFactors.set(outPos, lbFac.get(inPos));
                ++outPos;
            }
        }

        redFactors.set(numLMSSuff, true);

        //--------------------------------------------------------------------------------------------------------------
        // Encode the input string using labels for its LMS inf-suffixes to obtain the reduced version of the problem
        //--------------------------------------------------------------------------------------------------------------
        Tnum *redStr;

        try {
            redStr = new Tnum[numLMSSuff];
        }
        catch (const std::bad_alloc &e) {
            return -1;
        }

        for (Tnum inPos=len-1, outPos=numLMSSuff-1; inPos>=numLMSSuff; --inPos) {
            if (sa[inPos] !=0) {
                redStr[outPos] = sa[inPos] - 1;
                --outPos;
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        // Compute circular suffix array of the encoded string
        //--------------------------------------------------------------------------------------------------------------
        circularSuffixArray(redStr, sa, numLMSSuff, redFactors, numLabels);


        for (Tnum inPos=0, outPos=0; inPos<len; ++inPos) {
            if (isLMSPos(inPos, lbFac, suffType) && !spcSuff.get(inPos)) {
                redStr[outPos] = inPos;
                ++outPos;
            }
        }

        for (Tnum i=0; i<numLMSSuff; ++i) {
            sa[i] = redStr[sa[i]];
        }

        delete[] redStr;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Induce the result for the original problem
    //------------------------------------------------------------------------------------------------------------------

    memcpy(tmpBuckets, buckets, alphSize*sizeof(Tnum));

    for (Tnum i=numLMSSuff; i<len; ++i) {
        sa[i] = -1;
    }

    //------------------------------------------------------------------------------------------------------------------
    // Insert LMS ins-suffixes into proper buckets (starTdatag from bucket end)
    //------------------------------------------------------------------------------------------------------------------
    for (Tnum i=numLMSSuff-1; i>=0; --i) {
        Tnum j = sa[i];
        sa[i] = -1;
        sa[tmpBuckets[inStr[j]+1]-1] = j;
        --tmpBuckets[inStr[j]+1];
    }

    //---------------------------------------------------------------------------------------------
    // Insert L inf-suffixes into the proper bucket (starTdatag from the beginning of the bucket)
    //---------------------------------------------------------------------------------------------
    memcpy(tmpBuckets, buckets, alphSize*sizeof(Tnum));
    preSortSuffixexL(inStr, sa, len, lbFac, suffType, spcSuff, tmpBuckets, alphSize);

    //------------------------------------------------------------------------------------------------------------------
    // Insert S inf-suffixes into the proper bucket (starTdatag from the bucket end)
    //------------------------------------------------------------------------------------------------------------------
    memcpy(tmpBuckets, buckets, alphSize*sizeof(Tnum));
    preSortSuffixesS(inStr, sa, len, lbFac, suffType, tmpBuckets, alphSize);

    delete[] buckets;
    delete[] tmpBuckets;
    delete[] labelBuff;

    return 0;
}


#endif //_BBWT_INTERNAL_HPP_
