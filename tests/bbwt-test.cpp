/**
 * Bijective Burrows-Wheeler Transform computation testing.
 * Input data is read from a file, then BBWT and its inversion are computed
 * and compared to the input data.
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
#include <iomanip>
#include <cstdio>
#include <chrono>
#include <new>

#include "bbwt.hpp"

using namespace std;
using Tnum = int;


int main(int argc, char **argv) {
    unsigned char *inData = nullptr;
    unsigned char *bbwtData = nullptr;
    unsigned char *outData = nullptr;
    Tnum *csa = nullptr;
    FILE *inFile;

    if(argc != 2) {
        cerr << "Usage " << argv[0] << " input_file" << endl;

        return 1;
    }

    //-------------------------------------------------------------------------
    // Read data from the input file
    //-------------------------------------------------------------------------

    inFile = fopen(argv[1], "rb");

    fseek(inFile, 0, SEEK_END);
    Tnum dataSize = ftell(inFile);
    rewind(inFile);

    cout << "-- Input size = " << dataSize << " B --" << endl;

    try {
        inData = new unsigned char[dataSize];
        bbwtData = new unsigned char[dataSize];
        outData = new unsigned char[dataSize];
        csa = new Tnum[dataSize];
    }
    catch (const bad_alloc &e) {
        cerr << argv[0] << ": Memory allocation error" << endl;

        return 2;
    }

    Tnum dataCount = fread((char*) inData, sizeof(char), dataSize, inFile);
    fclose(inFile);

    if (dataCount != dataSize) {
        cerr << argv[0] << " error: input data read partially" << endl;

        return 1;
    }

    //-------------------------------------------------------------------------
    // Compute BBWT and measure the computation time
    //-------------------------------------------------------------------------

    cout << "-- Computing BBWT --" << endl;
    auto start = chrono::high_resolution_clock::now();

    if (bbwt(inData, bbwtData, csa, dataSize) != 0) {
        cerr << argv[0] << " error: BBWT computation failed" << endl;

        return -1;
    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << "-- Runtime " << duration.count()/1000 << "." << fixed << setprecision(3) << duration.count()%1000 << " s --" << endl;


    //-------------------------------------------------------------------------
    // Inverse the BBWT and compare the result to the input data
    //-------------------------------------------------------------------------

    cout << "-- Computing BBWT inversion --" << endl;

    if (unbbwt(bbwtData, outData,dataSize) != 0) {
        cerr << argv[0] << " error: inverse BBWT computation failed" << endl;

        return -1;
    }

    for (Tnum pos=0; pos < dataSize; ++pos) {
        if (inData[pos] != outData[pos]) {
            cout << "\t" << pos << ": [" << inData[pos] << "] != [" << outData[pos] << "]" << endl;
        }
    }

    cout << "-- Finished --" << endl;

    delete[] inData;
    delete[] bbwtData;
    delete[] outData;
    delete[] csa;

    //-------------------------------------------------------------------------

    return 0;
}
