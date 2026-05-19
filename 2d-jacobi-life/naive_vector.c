#include "define.h"
void naive_vector(int * A, int NX, int NY, int T){
	int t, x, y;
	int (* A_correct)[NX + 2 * XSTART][ NY + 2 * YSTART] =  (int (*)[NX + 2 * XSTART][ NY + 2 * YSTART]) A;
	for (t = 0; t < T; t++) {
		for (x = XSTART; x < NX + XSTART; x++) {
			#pragma vector always
			#pragma ivdep
			for (y = YSTART; y < NY + YSTART; y++) {
				Compute_scalar(A_correct, t, x, y);
			}
		}
	}
}
