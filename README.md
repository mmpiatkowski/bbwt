
## Linear-time Bijective Burrows-Wheeler Transform Construction

This repository contains the implementation of the linear-time Bijective Burrows-Wheeler Transform
construction algorithm described in:

[1] Hideo Bannai, Juha Kärkkäinen, Dominik Köppl, Marcin Piątkowski 
*Constructing the Bijective and the Extended Burrow-Wheeler Transform in Linear Time*,
in Proceeding of 32nd Annual Symposium on Combinatorial Pattern Matching (CPM 2021)

The implementation has been intended as a proof of concept rather than the most efficient one.
Therefore, we focused mainly on simplicity and ease of understanding over the speed of computation.


## Interface

The implementation provides the following functions:

* Linear-time Circular Suffix Array computation based ons SAIS induced sorting algorithm

```c++
/** Computes the circular suffix array of inStr.
 * @param inStr input data
 * @param csa buffer where computed circular suffix array is stored
 * @param len size of the input data
 * @param alphSize size of the alphabet
 * @return 0 after successful computation, non-zero in case of any error */
template<typename Tdata, typename Tnum>
int circularSuffixArray(const Tdata *inStr, Tnum *csa, Tnum len, const Tnum alphSize = 256);
```

* Linear-time Bijective Burrows-Wheeler construction 

```c++
/** Computes Bijective Burows-Wheeler Transform of inStr.
* @param inStr input data
* @param outStr buffer where the computed BBWT is stored
* @param len the size of the input data
* @param alphSize size of the alphabet
* @return 0 after successful computation, non-zero in case of any error */
template<typename Tdata, typename Tnum>
int bbwt(const Tdata *inStr, Tdata *outStr, Tnum *csa, Tnum len, const Tnum alphSize = 256);
```

* Inverse of Bijective Burrows-Wheeler Transform

```c++
/** Computes the inverse of Bijective Burrows-Wheeler Transform of inStr.
* @param inStr input data
* @param outStr buffer where the computed inverse of BBWT is stored
* @param len the size of the input data
* @param alphSize size of the alphabet
* @return 0 after successful computation, non-zero in case of any error */
template<typename Tdata, typename Tnum>
int unbbwt(const Tdata *inStr, Tdata *outStr, Tnum len, const Tnum alphSize = 256);
```

## Usage

The circular suffix array of a given text may be computed as follows:

```c++
#include "bbwt.hpp"
unsigned char * text = ...;
int * csa = ...; // Use long if the input size exceeds std::numeric_limits<int>::max()
int length = ...; // Size of the text
...
circularSuffixArray(text, csa, length); 
```

The Bijective Burrows-Wheeler Transform of a given text may be computed as follows:

```c++
#include "bbwt.hpp"
unsigned char * text = ...;
unsigned char * output = ...;
int * csa = ...; // Use long if the input size exceeds std::numeric_limits<int>::max()
int length = ...; // Size of the text
...
bbwt(text, output, csa, length); 
```

## Examples

We provided the following example programs:
* **bbwt-main.cpp** - Computation of BBWT for data read from a file.
  The result is stored in a file.
* **bbwt-console.cpp** - Computation of BBWT for the data read from the standard input.
  The result is printed to standard output.
* **csa-console.cpp** - Computation of circular suffix array for the data read from the standard input.
  The result is printed to the standard output.
  
  
## Experimental results

The abbreviations in the tables below are as follows:
* File - file name
* S - file size (in Bytes)
* A - alphabet size (number of distinct characters)
* \#LF - the number of Lyndon factors
* \#ULF - the number of unique Lyndon factos
* BBWT - the number of character runs in Bijective Burrows-Wheeler Transform
* BWT - the number of character runs in Burrows-Wheler Transform


### Calgary Corpus

| File   |      S |   A |  \#LF | \#ULF |   BBWT |    BWT |
|:-------|-------:|----:|------:|------:|-------:|-------:|
| bib    | 111261 |  81 |     6 |     6 |  36971 |  36964 |
| book1  | 768771 |  82 |    12 |    12 | 386264 | 386263 |
| book2  | 610856 |  96 |    27 |    27 | 239378 | 239367 |
| geo    | 102400 | 256 |    20 |     8 |  65781 |  65778 |
| news   | 377109 |  98 |    24 |    24 | 158607 | 158592 |
| obj1   |  21504 | 256 |   991 |     6 |  10616 |  10616 |
| obj2   | 246814 | 256 |    10 |    10 |  78814 |  78814 |
| paper1 |  53161 |  95 |     9 |     9 |  22146 |  22140 |
| paper2 |  82199 |  91 |    16 |    16 |  36689 |  36687 |
| paper3 |  46526 |  84 |    14 |    14 |  22569 |  22566 |
| paper4 |  13286 |  80 |     6 |     6 |   6904 |   6903 |
| paper5 |  11954 |  91 |     6 |     6 |   5938 |   5935 |
| paper6 |  38105 |  93 |    15 |    15 |  16048 |  16046 |
| pic    | 513216 | 159 | 36319 |     4 |  64691 |  64690 |
| progc  |  39611 |  92 |    12 |    12 |  15709 |  15707 |
| progl  |  71646 |  87 |    77 |     7 |  19446 |  19442 |
| progp  |  49379 |  89 |    12 |    12 |  12825 |  12823 |
| trans  |  93695 |  99 |   228 |    13 |  19456 |  19453 |


### Canterbury Corpus

| File         |       S |   A |  \#LF | \#ULF |   BBWT |    BWT |
|:-------------|--------:|----:|------:|------:|-------:|-------:|
| alice29.txt  |  152089 |  74 |     3 |     3 |  66903 |  66902 |
| asyoulik.txt |  125179 |  68 |     2 |     2 |  62366 |  62364 |
| cp.html      |   24603 |  86 |     8 |     8 |   9201 |   9198 |
| fields.c     |   11150 |  90 |    13 |    13 |   3417 |   3409 |
| grammar.lsp  |    3721 |  76 |     8 |     6 |   1340 |   1344 |
| kennedy.xls  | 1029744 | 256 |     9 |     9 | 234842 | 234838 |
| lcet10.txt   |  426754 |  84 |     6 |     6 | 165712 | 165709 |
| plrabn12.txt |  481861 |  81 |     6 |     6 | 243558 | 243557 |
| ptt5         |  513216 | 159 | 36319 |     4 |  64691 |  64690 |
| sum          |   38240 | 255 |    13 |    10 |  13262 |  13262 |
| xargs.1      |    4227 |  74 |     9 |     9 |   2009 |   2008 |
              

### Silesia Corpus

| File    |     Size |   A |  \#LF | \#ULF |     BBWT |      BWT |
|:--------|---------:|----:|------:|------:|---------:|---------:|
| dickens | 10192446 | 100 |    17 |    17 |  4374629 |  4374598 |
| mozilla | 51220480 | 256 |  8322 |    16 | 19498071 | 19498064 |
| mr      |  9970564 | 256 |   179 |    23 |  3444892 |  3444884 |
| nci     | 33553445 |  62 |     6 |     6 |  2195819 |  2195816 |
| ooffice |  6152192 | 256 | 15201 |     8 |  3215176 |  3215179 |
| osdb    | 10085684 | 256 |    15 |    15 |  3224389 |  3224386 |
| reymont |  6627202 | 256 |    26 |    26 |  1935978 |  1935966 |
| samba   | 21606400 | 256 | 10304 |    17 |  5220562 |  5220551 |
| sao     |  7251944 | 256 |     6 |     6 |  5524741 |  5524740 |
| x-ray   |  8474240 | 256 |     7 |     7 |  4796213 |  4796213 |
| xml     |  5345280 | 104 |  5955 |     8 |   581963 |   581955 |


### Pizza \& Chili Corpus

| File     |          S |   A | \#LF | \#ULF |      BBWT |        BWT |
|:---------|-----------:|----:|-----:|------:|----------:|-----------:|
| dblp.xml |  296135874 |  97 |   15 |    15 |  41037558 |   41037553 |
| dna      |  403927746 |  16 |   18 |    18 | 243492872 |  243492866 |
| english  | 2210395553 | 239 |   18 |    18 | 658301004 |  658301008 |
| pitches  |   55832855 | 133 |   39 |    39 |  23430040 |   23429976 |
| proteins | 1184051855 |  27 |   30 |    30 | 441858493 |  441858468 |
| sources  |  210866607 | 230 |   31 |    31 |  47896880 |   47896806 |


### Pizza \& Chili Repetitive Corpus

| File             |         S |   A | \#LF | \#ULF |     BBWT |      BWT |
|:-----------------|----------:|----:|-----:|------:|---------:|---------:|
| fib41            | 267914296 |   2 |   21 |    21 |       41 |        3 |
| rs.13            | 216747218 |   2 |   27 |    27 |      123 |       75 |
| tm29             | 268435456 |   2 |   41 |    41 |       81 |       81 |
| dblp.xml.00001.1 | 104857600 |  89 |    7 |     7 |   172500 |   172487 |
| dblp.xml.00001.2 | 104857600 |  89 |   23 |    23 |   175626 |   175616 |
| dblp.xml.0001.1  | 104857600 |  89 |    9 |     9 |   240550 |   240533 |
| dblp.xml.0001.2  | 104857600 |  89 |   21 |    21 |   270213 |   270203 |
| dna.001.1        | 104857600 |   5 |   18 |    18 |  1716857 |  1716806 |
| english.001.2    | 104857600 | 106 |   29 |    29 |  1449562 |  1449517 |
| proteins.001.1   | 104857600 |  21 |   19 |    19 |  1278237 |  1278199 |
| sources.001.2    | 104857600 |  98 |   50 |    50 |  1213519 |  1213426 |
| Escherichia_Coli | 112689515 |  15 |   13 |    13 | 15044536 | 15044485 |
| cere             | 461286644 |   5 |   21 |    21 | 11574705 | 11574639 |
| coreutils        | 205281778 | 236 |   17 |    17 |  4684513 |  4684458 |
| einstein.de.txt  |  92758441 | 117 |   21 |    21 |   101391 |   101369 |
| einstein.en.txt  | 467626544 | 139 |   59 |    59 |   290279 |   290237 |
| influenza        | 154808555 |  15 |   10 |    10 |  3022821 |  3022820 |
| kernel           | 257961616 | 160 |   32 |    32 |  2791456 |  2791366 |
| para             | 429265758 |   5 | 1238 |    22 | 15636838 | 15636738 |
| world_leaders    |  46968181 |  89 |   13 |    12 |   573506 |   573485 |

