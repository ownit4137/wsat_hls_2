#include "wsat_hls.h"
#define CSIME

#ifndef CSIM
#include "hls_print.h"
#endif

extern "C" {

inline int psrandom(int& seed) {
	ap_uint<32> lfsr = seed;
	bool b_32 = lfsr.get_bit(32-32);
	bool b_22 = lfsr.get_bit(32-22);
	bool b_2 = lfsr.get_bit(32-2);
	bool b_1 = lfsr.get_bit(32-1);
	bool new_bit = b_32 ^ b_22 ^ b_2 ^ b_1;
	lfsr = lfsr >> 1;
	lfsr.set_bit(31, 0);
	lfsr.set_bit(30, new_bit);
	seed = lfsr.to_int();
	return lfsr.to_int();
}

inline uint8_t psrandom8b(int& seed) {
	uint8_t r = psrandom(seed);
	return r;
}

struct reqctrl {
	bool isloc;
	bool isinit;
	int id;
	int length;
};

struct rspctrl {
	int val[DSIZE];
	int size;
	int numcls;
};

struct changes {
	int id;
	bool plus;
};

#ifdef CSIM

// int ClausesVec_c[MAXNCLS * K];
hls::vector<int, DSIZE> ClausesVec_c[MAXNCLS * K / DSIZE];
// int VarsLocVec_c[MAXNLIT * DSIZE * R];
hls::vector<int, DSIZE> VarsLocVec_c[MAXNLIT * R];


void forcsim(hls::vector<int, DSIZE> ClausesVec[], hls::vector<int, DSIZE> VarsLocVec[]) {
	for (int i = 0; i < MAXNLIT * R; i++) {
		for (int j = 0; j < DSIZE; j++) {
			VarsLocVec_c[i][j] = VarsLocVec[i][j];
		}
	}

	for (int i = 0; i < MAXNCLS * K / DSIZE; i++) {
		for (int j = 0; j < DSIZE; j++) {
			ClausesVec_c[i][j] = ClausesVec[i][j];
		}
	}
}

bool verify(bool VATArr[], int numClauses) {
	bool result = true;
	verify_loop: for (int c = 1; c <= numClauses; c++) {
		bool cls = false;

		for (int k = 0; k < K; k++) {
			for (int i = 0; i < DSIZE; i++) {
				int lit = ClausesVec_c[c * K / DSIZE + k / DSIZE][i];
				if (lit == 0) break;
				if (VATArr[ABS(lit)] == (lit > 0)) cls = true;
			}
		}
		if (cls == false) {
			result = false;
			break;
		}
	}
	return result;
}

void YalSATread3(hls::stream<reqctrl>& reqStream, hls::stream<rspctrl>& rspStream) {
	reqctrl r;
	rspctrl w;

	int tempcls[K];
	r = reqStream.read();

	int var = r.id;

	if ((var == -1) && (r.isloc == false)) {
		return;
	}

	if (r.isloc) {
		int numD = r.length;				// 1~4: 1, 5~8: 2

		for (int d = 0; d < numD; d++) {	// one line
			for (int i = 0; i < DSIZE; i++) {
				// int idx = GETPOS(var) * R * DSIZE + d * DSIZE + i;
				int idx = GETPOS(var) * R + d;
				int loc = VarsLocVec_c[idx][i];
				w.val[i] = loc;
			}

			w.size = DSIZE;
			rspStream.write(w);
		}

		// std::cout << "read for var " << var << " " << numD << " blocks are sent\n";
	}
	else if (r.isinit) {
		initloop: for (int c = 1; c <= r.length; c++) {
			int cnt = 0;

			w.numcls = c;
			cls_io: for (int k = 0; k < K / DSIZE; k++) {
				hls::vector<int, DSIZE> t = ClausesVec_c[c * K / DSIZE + k];
				int cnt = 0;
				for (; cnt < DSIZE; cnt++) {
					int cls = t[cnt];
					w.val[cnt] = t[cnt];
					if (cls == 0) break;
				}
				w.size = cnt;

				if (cnt != 0)	rspStream.write(w);
			}
		}
	}
	else {

		int endK = r.length;
		int numD = ((endK - 1) / DSIZE) + 1;				// 1~4: 1, 5~8: 2

		for (int d = 0; d < numD; d++) {
			w.size = (d == numD - 1) ? ((endK - 1) % DSIZE) + 1 : DSIZE;						// mod 2
			for (int i = 0; i < w.size; i++) {
				w.val[i] = ClausesVec_c[var * K / DSIZE + d][i];
			}
			rspStream.write(w);
			// std::cout << "for " << c << " cnt " << cnt << " block " << d << " size " << w.size << "\n";
		}
	}
	// std::cout << "Read finished " << var << "\n";
}
#endif


void YalSATread2(hls::stream<reqctrl>& reqStream, hls::stream<rspctrl>& rspStream,
		hls::vector<int, DSIZE> VarsLocVec[], hls::vector<int, DSIZE> ClausesVec[]) {

	std::cout << "read starts" << "\n";
	// hls::print("read starts\n");

	reqctrl r;
	r = reqStream.read();
	if (not r.isinit) {return;}

	rspctrl w;

	initloop: for (int c = 1; c <= r.length; c++) {
#pragma HLS PIPELINE II = 1
		int cnt = 0;

		w.numcls = c;
		cls_io: for (int k = 0; k < K / DSIZE; k++) {
			hls::vector<int, DSIZE> t = ClausesVec[c * K / DSIZE + k];
			int cnt = 0;
			for (; cnt < DSIZE; cnt++) {
				int cls = t[cnt];
				w.val[cnt] = cls;
				if (cls == 0) break;
			}
			w.size = cnt;
			if (cnt != 0) rspStream.write(w);
		}
	}

	ioloop: while (true) {
		r = reqStream.read();
		int var = r.id;
		if ((var == -1) && (r.isloc == false)) {
			break;
		}

		if (r.isloc) {
			int numD = r.length;				// 1~4: 1, 5~8: 2

			loc_blocks: for (int d = 0; d < numD; d++) {	// one line
#pragma HLS PIPELINE II = 1

				hls::vector<int, DSIZE> t = VarsLocVec[GETPOS(var) * R + d];
				loc_dsize: for (int i = 0; i < DSIZE; i++) {
					w.val[i] = t[i];
				}

				w.size = DSIZE;
				rspStream.write(w);
			}

		}
		else {				// uc
			int endK = r.length;
			int numD = ((endK - 1) / DSIZE) + 1;				// 1~4: 1, 5~8: 2

			uc_blocks: for (int d = 0; d < numD; d++) {
#pragma HLS PIPELINE II = 1
				hls::vector<int, DSIZE> t = ClausesVec[var * K / DSIZE + d];
				w.size = (d == numD - 1) ? ((endK - 1) % DSIZE) + 1 : DSIZE;						// mod 2

				for (int i = 0; i < w.size; i++) {
					w.val[i] = t[i];
				}
				rspStream.write(w);
				// std::cout << "for " << c << " cnt " << cnt << " block " << d << " size " << w.size << "\n";
			}
		}
	}
}


void YalSATmain(hls::stream<reqctrl>& reqStream, hls::stream<rspctrl>& rspStream,
		hls::vector<short, SDSIZE> varND_off[], hls::vector<short, SDSIZE> clsK_off[],
		int numVars, int numClauses, int maxflip, bool issolved, int s) {
	// std::cout << "main starts" << "\n";
#ifdef CSIM
	clock_t start = clock();
#endif

	char ClausesCost2[MAXNCLS / DSIZE + 1][DSIZE];
	bool VATArr[MAXNVAR];
	int UCBArr[UCBSIZE];

	int probsLT[160] = {40710, 8734, 3655, 1985, 1240, 846, 612, 463, 362, 291, 238, 199, 168, 144, 125, 109, 96, 86, 76, 69,
			62, 56, 51, 47, 43, 40, 37, 34, 32, 29, 27, 26, 24, 23, 21, 20, 19, 18, 17, 16,
			15, 14, 14, 13, 12, 12, 11, 11, 10, 10, 9, 9, 9, 8, 8, 8, 7, 7, 7, 7,
			6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2,
			2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
	short numcrit[MAXNVAR];

	int XOR_UCBidx2[MAXNCLS / DSIZE + 1][DSIZE];
	int UCBlast;
	reqctrl req;
	rspctrl rsp;

	short varND[MAXNLIT];
	short clsK[MAXNCLS];
	int seed = s;
	int NC = numClauses;
	int NV = numVars;

	changes UCBlist[DSIZE][R * 2];
	short UCBlistidx[DSIZE];
	changes numcritlist[DSIZE][R * 2];
	short numcritlistidx[DSIZE];

	init_idx: for (int i = 0; i < DSIZE; i++) {
#pragma HLS PIPELINE II = 1
		UCBlistidx[i] = 0;
		numcritlistidx[i] = 0;
	}

	copy_nd: for (int i = 0; i <= NV * 2 / SDSIZE + 1; i++) {
#pragma HLS PIPELINE II = 1

		hls::vector<short, SDSIZE> t = varND_off[i];
		for (int j = 0; j < SDSIZE; j++) {
			varND[i * SDSIZE + j] = t[j];
		}
	}

	copy_k: for (int i = 0; i <= NC / SDSIZE + 1; i++) {
#pragma HLS PIPELINE II = 1

		hls::vector<short, SDSIZE> t = clsK_off[i];
		for (int j = 0; j < SDSIZE; j++) {
			clsK[i * SDSIZE + j] = t[j];
		}
	}

#pragma HLS ARRAY_PARTITION variable=XOR_UCBidx2 complete dim=2
#pragma HLS ARRAY_PARTITION variable=ClausesCost2 complete dim=2
#pragma HLS ARRAY_PARTITION variable=numcritlist complete dim=1
#pragma HLS ARRAY_PARTITION variable=UCBlist complete dim=1

	var: for (int i = 0; i <= NV; i++) {
#pragma HLS PIPELINE II = 1
		VATArr[i] = psrandom(seed) % 2 == 0 ? true : false;		// on-chip
		// VATArr[i] = i % 2 == 0 ? true : false;		// on-chip
		numcrit[i] = 0;
	}
	UCBlast = 0;

	req.isloc = false;
	req.isinit = true;
	req.length = NC;
	reqStream.write(req);

#ifdef CSIM
	YalSATread3(reqStream, rspStream);
#endif

	cls: for (int c = 1; c <= NC; c++) {
#pragma HLS PIPELINE II = 1

		int row = c / DSIZE;
		int col = c % DSIZE;

		XOR_UCBidx2[row][col] = 0;

		int endK = clsK[c];
		int numD = ((endK - 1) / DSIZE) + 1;

		int cost = 0;
		int tl;
		// std::cout << "C" << c << " " << endK << " [";

		clsinit_blocks: for (int d = 0; d < numD; d++) {
#pragma HLS PIPELINE II = 1
			rsp = rspStream.read();
			int size = rsp.size;


			clsinit_dsize: for (int i = 0; i < size; i++) {
				int lit = rsp.val[i];
				if (VATArr[ABS(lit)] == (lit > 0)) {
					cost++;
					tl = lit;
					XOR_UCBidx2[row][col] ^= ABS(lit);
					// std::cout << lit << " ";
				}
			}
		}



		if (cost == 0) {
			UCBArr[UCBlast] = c;			// on-chip: int
			XOR_UCBidx2[row][col] = UCBlast++;		// on-chip: int
			// std::cout << " cost 0" << UCBlast << " \n";
		}
		else if (cost == 1) {
			numcrit[ABS(tl)]++;
			// std::cout << " cost 1 nc" << ABS(tl) << ":" << numcrit[ABS(tl)] << "\n";
		}
		else {
			// std::cout << " cost " << cost << "\n";
		}
		ClausesCost2[row][col] = cost;
	}

/////////////////////// [[  FLIP  ]] //////////////////////////

	flip: for (unsigned long long f = 0; f < MAX_FLIPS; f++) {
#pragma HLS PIPELINE off


		if (f % (MAX_FLIPS / PRINT_FLIP_TIMES) == 0) {
#ifdef CSIM
			std::cout << "flip-" << f << " " << UCBlast << " remain \n";
#endif
		}


		/////////////////////// 1. IO uc //////////////////////////

		int ucnum = UCBArr[psrandom(seed) % UCBlast];
		int var_flip;

		req.isloc = false;
		req.isinit = false;
		req.id = ucnum;
		req.length = clsK[ucnum];
		reqStream.write(req);

#ifdef CSIM
		YalSATread3(reqStream, rspStream);
#endif

		// std::cout << "ucnum:" << ucnum << " endK:" << req.length << " ";

		int tempcls[K];
		int probs[K];
		int sumProb = 0;
		int endK = req.length;
		int numD = ((endK - 1) / DSIZE) + 1;				// 1~4: 1, 5~8: 2

		ucb_read: for (int d = 0; d < numD; d++) {
#pragma HLS PIPELINE II = 1
			rsp = rspStream.read();
			ucb_read_dsize: for (int i = 0; i < rsp.size; i++) {
				tempcls[d * DSIZE + i] = rsp.val[i];
			}
		}

		lookup_break: for (int i = 0; i < endK; i++) {
#pragma HLS PIPELINE II = 1
#pragma HLS loop_tripcount max=3

			short bv = numcrit[ABS(tempcls[i])];
			sumProb += probsLT[bv];
			probs[i] = sumProb;
			// std::cout << tempcls[i] << "_" << probs[i] << "(" << bv << ") ";
		}

		int randPosition = psrandom8b(seed) * sumProb / 256;

		// std::cout << randPosition << "r ";
		choose_var: for (int i = 0; i < endK; i++) {
#pragma HLS PIPELINE II = 1
#pragma HLS loop_tripcount max=3

			if (probs[i] >= randPosition) {
				var_flip = tempcls[i];
				// std::cout << "   " << i << " " << var_flip << " chosen\n";
				break;
			}
		}



		/////////////////////// 2. IO loc //////////////////////////

		req.isloc = true;
		req.id = var_flip;
		int d1 = varND[GETPOS(var_flip)];
		req.length = d1;
		reqStream.write(req);

#ifdef CSIM
		YalSATread3(reqStream, rspStream);
#endif

		req.isloc = true;
		req.id = -var_flip;
		int d2 = varND[GETPOS(-var_flip)];
		req.length = d2;
		reqStream.write(req);

#ifdef CSIM
		YalSATread3(reqStream, rspStream);
#endif
		// std::cout << var_flip << "chosen\n   d1[" ;

		/////////////////////// 3. save //////////////////////////


		tl_inc_rsp_loop: for (int t = 0; t < d1; t++) {
#pragma HLS dependence variable=ClausesCost2 type=inter false
#pragma HLS dependence variable=XOR_UCBidx2 type=inter false
#pragma HLS PIPELINE II = 1



			rsp = rspStream.read();
			one_line1: for (int i = 0; i < DSIZE; i++) {
#pragma HLS UNROLL
				int cn = rsp.val[i];
				// std::cout << cn << " ";

				if (cn > 0) {

					int row = cn / DSIZE;
					int cost = ClausesCost2[row][i];
					int critv;

					if (cost == 0) {
						critv = XOR_UCBidx2[row][i];

						changes p;
						p.id = cn;
						p.plus = false;
						UCBlist[i][UCBlistidx[i]++] = p;		// 0 -> 1 ucberase;

						p.id = ABS(var_flip);
						p.plus = true;
						numcritlist[i][numcritlistidx[i]++] = p;

					} else if (cost == 1) {
						critv = XOR_UCBidx2[row][i];

						changes p;
						p.id = critv;
						p.plus = false;
						numcritlist[i][numcritlistidx[i]++] = p;
						critv = critv ^ ABS(var_flip);
					} else {
						critv = XOR_UCBidx2[row][i] ^ ABS(var_flip);
					}

					XOR_UCBidx2[row][i] = critv;
					ClausesCost2[row][i] = cost + 1;		// can be moved  c+1
				}
			}
			// std::cout << "\n";
		}

		// std::cout << "] \n   d2[";
		tl_dec_rsp_loop: for (int t = 0; t < d2; t++) {
#pragma HLS dependence variable=ClausesCost2 type=inter false
#pragma HLS dependence variable=XOR_UCBidx2 type=inter false
#pragma HLS PIPELINE II = 1

			rsp = rspStream.read();
			one_line2: for (int i = 0; i < DSIZE; i++) {
#pragma HLS UNROLL
				int cn = rsp.val[i];
				// std::cout << cn << " ";

				if (cn > 0) {
					int row = cn / DSIZE;
					int critv = XOR_UCBidx2[row][i];
					int cost = ClausesCost2[row][i];
					XOR_UCBidx2[row][i] = critv ^ ABS(var_flip);

					if (cost == 1) {
						changes p;
						p.id = cn;
						p.plus = true;
						UCBlist[i][UCBlistidx[i]++] = p;		// 1-> 0 ucbinsert

						p.id = ABS(var_flip);
						p.plus = false;
						numcritlist[i][numcritlistidx[i]++] = p;

					} else if (cost == 2) {
						changes p;
						p.id = critv ^ ABS(var_flip);
						p.plus = true;
						numcritlist[i][numcritlistidx[i]++] = p;
					}

					ClausesCost2[row][i] = cost - 1;		// can be moved
				}
			}
			// std::cout << "\n";
		}

		// std::cout << "] \n";

		/////////////////////// 4. update //////////////////////////


		fifo_update1_2: for (int i = 0; i < DSIZE; i++) {
			numcrit_update: for (int j = 0; j < numcritlistidx[i]; j++) {
#pragma HLS PIPELINE II = 1
				changes crit = numcritlist[i][j];
				if (crit.id > 0) {
					if (crit.plus) {
						numcrit[crit.id]++;
					} else {
						numcrit[crit.id]--;
					}
				}
			}
			numcritlistidx[i] = 0;
		}

		fifo_update2: for (int i = 0; i < DSIZE; i++) {
#pragma HLS PIPELINE II = 1
			ucb_update: for (int j = 0; j < UCBlistidx[i]; j++) {
				changes cand = UCBlist[i][j];
				int cn = cand.id;
				if (cn > 0) {
					if (cand.plus) {	// insert
						UCBArr[UCBlast] = cn;
						XOR_UCBidx2[cn / DSIZE][cn % DSIZE] = UCBlast++;
					}
					else {		// erase
						int lastElem = UCBArr[--UCBlast];

						int outidx = XOR_UCBidx2[cn / DSIZE][cn % DSIZE];
						UCBArr[outidx] = lastElem;
						XOR_UCBidx2[lastElem / DSIZE][lastElem % DSIZE] = outidx;
						XOR_UCBidx2[cn / DSIZE][cn % DSIZE] = ABS(var_flip);
					}
				}
			}
			UCBlistidx[i] = 0;
		}

		VATArr[ABS(var_flip)] = 1 - VATArr[ABS(var_flip)];

		// finish cond
		if (UCBlast == 0){
			issolved = true;
			maxflip = f;
			break;
		}
	}

	// end flip

	req.isloc = false;
	req.isinit = false;
	req.id = -1;

	reqStream.write(req);

#ifdef CSIM
	// verify answer
	YalSATread3(reqStream, rspStream);
	// std::cout << reqStream.size() << " " << rspStream.size() << "\n";

	if (issolved) {
		bool v = verify(VATArr, numClauses);
		if (v) { std::cout << "Solver found a solution. Verified. | "; }
		else {
			issolved = false;
			std::cout << "Wrong solution. UCB count : " << UCBlast << "\n";
		}
	}
	else {
		std::cout << "Solver could not find a solution. | ";
		maxflip = MAX_FLIPS;
	}
#endif

#ifdef CSIM
	clock_t end = clock();
	std::cout << "Solver completed in: " << (double)(end - start)/CLOCKS_PER_SEC << " seconds | flip: " << maxflip << std::endl;
#else
	std::cout << "Solver completed in: " << 0 << " seconds | flip: " << maxflip << std::endl;
#endif
}


void yalsat(hls::vector<short, DSIZE * 2> varND_off[], hls::vector<short, DSIZE * 2> clsK_off[],
		hls::vector<int, DSIZE> ClausesVec[], hls::vector<int, DSIZE> VarsLocVec[],
		int numVars, int numClauses, int maxflip, bool issolved, int s) {
	//std::cout << "yalsat start\n";

	hls::stream<reqctrl> reqStream;
	hls::stream<rspctrl> rspStream;

#ifndef __SYNTHESIS__
	//std::cout << "nv: " << numVars << " nc: " << numClauses << "\n";
#ifdef CSIM
	forcsim(ClausesVec, VarsLocVec);
	YalSATmain(reqStream, rspStream, varND_off, clsK_off, numVars, numClauses, maxflip, issolved, s);
	return;
#endif
#endif

//	hls::vector<int, DSIZE> ClausesVec[MAXNCLS * K / DSIZE];
//	hls::vector<int, DSIZE> VarsLocVec[MAXNLIT * R];
//	hls::vector<short, SDSIZE> clsK_v[MAXNCLS / SDSIZE + 1];
//	hls::vector<short, SDSIZE> varND_v[MAXNLIT / SDSIZE + 1];





#pragma HLS INTERFACE mode=m_axi bundle=m0 port=varND_off depth = 101
#pragma HLS INTERFACE mode=m_axi bundle=m1 port=clsK_off depth = 101
#pragma HLS INTERFACE mode=m_axi bundle=m2 port=ClausesVec depth = 6400
#pragma HLS INTERFACE mode=m_axi bundle=m3 port=VarsLocVec depth = 512000


#pragma HLS DATAFLOW
#pragma HLS STREAM variable=reqStream depth=32
#pragma HLS STREAM variable=rspStream depth=32

	YalSATmain(reqStream, rspStream, varND_off, clsK_off, numVars, numClauses, maxflip, issolved, s);
	YalSATread2(reqStream, rspStream, VarsLocVec, ClausesVec);
}

}
