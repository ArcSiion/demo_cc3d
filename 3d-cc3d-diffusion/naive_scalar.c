#include "define.h"

void naive_scalar(float *A, int NX, int NY, int NZ, int T) {
	float (*B)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART] =
		(float (*)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART])A;

	for (int t = 0; t < T; t++) {
		for (int x = XSTART; x < NX + XSTART; x++) {
			for (int y = YSTART; y < NY + YSTART; y++) {
				#pragma novector
				for (int z = ZSTART; z < NZ + ZSTART; z++) {
					Compute_scalar(B, t, x, y, z);
				}
			}
		}
	}
}
