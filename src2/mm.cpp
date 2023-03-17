#include "hls_vector.h"
#include "hls_stream.h"
#include "ap_int.h"
#include "mm.h"

const int WIDTH = 64/sizeof(DTYPE);

extern "C" {

void mm(DTYPE *A,  DTYPE *B, DTYPE *AB,   int N )
{
#pragma HLS INTERFACE mode=m_axi bundle=m0 port=A 
#pragma HLS INTERFACE mode=m_axi bundle=m1 port=B 
#pragma HLS INTERFACE mode=m_axi bundle=m1 port=AB 


	DTYPE AB_block[M][M];
#pragma HLS ARRAY_PARTITION dim=2 type=complete variable=AB_block

	DTYPE Bj[M];
#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=Bj

	ib_loop: for(int ib = 0; ib < N/M; ib++) {
		jb_loop: for(int jb = 0; jb < N/M; jb++) {
			init_i_loop: for(int i = 0; i < M; i++) {
#pragma HLS pipeline II = 1
				init_j_loop: for(int j = 0; j < M; j++) {
#pragma HLS unroll
					AB_block[i][j] = 0;
				}
			}

			kb_loop: for(int kb = 0; kb < N/M; kb++) {
				k_loop: for(int k = 0; k < M; k++) {
					readB_j_loop: for(int j = 0; j < M; j++) {
#pragma HLS pipeline II = 1
						DTYPE B_temp = B[(kb*M+k)*N+jb*M+j];
						Bj[j] = B_temp;
					}

					i_loop: for(int i = 0; i < M; i++) {
#pragma HLS pipeline II=1
						DTYPE Ai =  A[((ib*M+i)*N+kb*M)+k];
						j_loop: for(int j = 0; j < M; j++) {
#pragma HLS unroll
							AB_block[i][j] += Ai * Bj[j];
						}
					}
				}
			}

			writeAB_i_loop: for(int i = 0; i < M; i++) {
				writeAB_j_loop: for(int j = 0; j < M; j++) {
#pragma HLS pipeline II = 1
					AB[(ib*M+i)*N+jb*M+j] = AB_block[i][j];
				}
			}
		}
	}
}

}
