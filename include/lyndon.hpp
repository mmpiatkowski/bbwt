#ifndef _LYNDON_HPP_
#define _LYNDON_HPP_

/**
 * Lyndon factorisation implementation.
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


/** Lyndon factorisation based on Duval's algorithm.
 *
 * @param inStr Input data.
 * @param factors BitVector os size (length + 1).
 *        Starting position of each Lyndon factor of inStr is marked with bit value 1.
 *        Additionally, the position after the last character of inStr is also marked
 *        with bit value 1.
 * @param unique Similar to <i>factors</b>, but for each Lyndon factor only the starting
 *        position of its first occurrence is marked.
 * @param length The size of the input data (the length of inStr).
 *
 * @return The number of all Lyndon factors of inStr.
*/
template<typename Tdata, typename Tnum>
Tnum lyndonFactors(const Tdata *inStr, Tnum length, BitVector<Tnum> *factors = nullptr, BitVector<Tnum> *unique = nullptr) {
    Tnum i = 0; 
    Tnum numFactors = 0;

    while (i < length) {
        Tnum j = i + 1, k = i;

        // Compute the length of the longest segment consisting of repeated occurrences
        // of the same Lyndon word.

        while (j < length && inStr[k] <= inStr[j]) {
            if (inStr[k] < inStr[j])
                k = i;
            else
                k++;
            j++;
        }

        // For each factor store its starting position (factors).
        // Moreover, store the starting position of the first occurrence of each factor (unique).

        if (unique)
            unique->set(i, true);

        while(i <= k) {
            if (factors)
                factors->set(i, true);
            ++numFactors;
            i += j - k;
        }
    }

    if (factors)
        factors->set(length, true);

    return numFactors;
}


#endif //_LYNDON_HPP_
