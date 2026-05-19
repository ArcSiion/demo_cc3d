#include "define.h"


void naive_scalar(double* A, int NX, int NY, int NZ, int T) {
	double (* A_correct)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;
    int t, x, y, z;
    for (t = 0; t < T; t++) {
		for (x = XSTART; x < NX + XSTART; x++) {
			for (y = YSTART; y < NY + YSTART; y++) {
				#pragma novector
				for (z = ZSTART; z < NZ + ZSTART; z++) {
					Compute_scalar(A_correct,t,x,y,z);
				}
			}
		}
	}

}