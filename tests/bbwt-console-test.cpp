/**
 * Bijective Burrows-Wheeler computation example.
 * Input data is read from the standard input,
 * then BBWT and its inversion are computed.
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

#include <cstdio>
#include <cstring>

#include "bbwt.hpp"

using namespace std;
using Tnum = int;


const int MaxDataSize = 2048;


int main() {
    unsigned char inData[MaxDataSize], bbwtData[MaxDataSize], outData[MaxDataSize];
    Tnum csa[MaxDataSize];

    printf("INPUT > ");

    while (fgets((char *)inData, MaxDataSize, stdin)) {
        Tnum len = strlen((char *)inData) - 1;
        inData[len] = bbwtData[len] = outData[len] = '\0';

        bbwt(inData, bbwtData, csa, len);
        printf("BBWT  > %s\n", bbwtData);
        unbbwt(bbwtData, outData, len);
        printf("UNBBWT> %s\n", outData);

        printf("\nINPUT > ");
    }

    return 0;
}
