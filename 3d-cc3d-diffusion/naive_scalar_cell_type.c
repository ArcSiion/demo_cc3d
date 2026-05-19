#include "define.h"

void naive_scalar_cell_type(float *A, cc3d_cell_type_t *cellType, const float *diffCoef, const float *decayCoef,
							int NX, int NY, int NZ, int T) {
	float (*B)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART] =
		(float (*)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART])A;
	cc3d_cell_type_t (*cellTypeField)[NY + 2 * YSTART][NZ + 2 * ZSTART] =
		(cc3d_cell_type_t (*)[NY + 2 * YSTART][NZ + 2 * ZSTART])cellType;

	for (int t = 0; t < T; t++) {
		for (int x = XSTART; x < NX + XSTART; x++) {
			for (int y = YSTART; y < NY + YSTART; y++) {
				#pragma novector
				for (int z = ZSTART; z < NZ + ZSTART; z++) {
					Compute_scalar_cell_type(B, cellTypeField, diffCoef, decayCoef, t, x, y, z);
				}
			}
		}
	}
}
