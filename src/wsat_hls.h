#ifndef WSAT_HLS_H
#define WSAT_HLS_H

#include <iostream>
#include <fstream>
#include <string>
#include "ap_int.h"
#include "hls_vector.h"
#include "hls_stream.h"
//#include <math.h>
//#include <ctype.h>
//#include <chrono>
//#include <set>
//#include <ctime>

extern "C" {

#define MAX_TRIES 1
#define MAX_FLIPS 500000
#define UCBSIZE 50000
#define PRINT_FLIP_TIMES 1
#define K 32 // multiple of 16
#define R 160
#define MAXNCLS 3200
#define MAXNVAR 1600
#define MAXNLIT 2 * MAXNVAR
#define GETPOS(L) (2*ABS(L)-(L<0)-1)	// -1(0) 1(1) -2(2) 2(3) -3(4) 3(5)
#define ABS(a) (((a) < 0) ? (-a) : (a))
#define DSIZE 16
#define SDSIZE 32
#define DUP 16


// void yalsat(hls::vector<short, SDSIZE>* varND_off, hls::vector<short, SDSIZE>* clsK_off,
// 		hls::vector<int, DSIZE>* ClausesVec, hls::vector<int, DSIZE>* VarsLocVec,
// 		int numVars, int numClauses, int& maxflip, bool& issolved, int seed);

}

#endif
