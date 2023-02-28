
#ifndef _BIT_VECTOR_HPP_
#define _BIT_VECTOR_HPP_

/**
 * Bit vector implementation.
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

#include <iostream>
#include <new>
#include <cstring>


template<typename Tnum>
class BitVector {
public:
    explicit BitVector(Tnum size) : numBits(size) {
        data = new unsigned char[(numBits >> 3) + 1]();
    }

    virtual ~BitVector() {
        delete[] data;
    }

    void clear() {
        memset(data, 0, ((numBits >> 3) + 1) * sizeof(unsigned char));
    }

    inline Tnum size() const {
        return numBits;
    }

    inline bool get(Tnum pos) const {
        return !(pos < 0 || pos >= numBits) && (data[pos >> 3] >> (pos & 7)) & 1;
    }

    inline void set(Tnum pos, bool val) {
        if (pos < 0 || pos >= numBits)
            return;

        if(val) {
            data[pos >> 3] |= (unsigned char)(1U << (pos & 7));
        }
        else {
            data[pos >> 3] &= (unsigned char)(~(1U << (pos & 7)));
        }
    }

    Tnum next(Tnum pos) const {
        Tnum j;
        unsigned int c;

        pos += 1;
        j = pos >> 3;
        c = data[j] >> (pos & 7);

        if(c == 0) {
            pos += 8 - (pos & 7);
            j += 1;

            for(; (c = data[j]) == 0; ++j, pos += 8) { }
        }

        for(; (c & 1) == 0; ++pos, c >>= 1) { }

        return pos;
    }

    Tnum prev(Tnum pos) const {
        Tnum j;
        unsigned int c;

        pos -= 1;

        if(0 <= pos) {
            j = pos >> 3;
            c = (data[j] << (7 - (pos & 7))) & 0xff;

            if(c == 0) {
                pos -= (pos & 7) + 1;
                j -= 1;

                for(; (0 <= j) && ((c = data[j]) == 0); --j, pos -= 8) { }

                if(c == 0) { c = 128; }
            }

            for(; (c & 128) == 0; --pos, c <<= 1) { }
        }

        return pos;
    }

private:
    Tnum numBits;
    unsigned char *data;
};


#endif //_BIT_VECTOR_HPP_
