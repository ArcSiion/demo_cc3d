#include "../common.h"

#define XSTART 1
#define YSTART 1
#define ZSTART 1

#define XSLOPE 1
#define YSLOPE 1
#define ZSLOPE 1

#define CC3D_DIFFUSION_COEF 0.1f
#define CC3D_DECAY_COEF 0.0f
#define CC3D_NEIGHBOR_COEF CC3D_DIFFUSION_COEF
#define CC3D_CENTER_COEF (1.0f - CC3D_DECAY_COEF - 6.0f * CC3D_DIFFUSION_COEF)

#define INIT ((float)(rand() % 1024))

#define Compute_scalar(A, t, x, y, z) \
	(A)[((t) + 1) % 2][x][y][z] = \
		CC3D_CENTER_COEF * (A)[(t) % 2][x][y][z] + \
		CC3D_NEIGHBOR_COEF * ( \
			(A)[(t) % 2][(x) - 1][y][z] + \
			(A)[(t) % 2][(x) + 1][y][z] + \
			(A)[(t) % 2][x][(y) - 1][z] + \
			(A)[(t) % 2][x][(y) + 1][z] + \
			(A)[(t) % 2][x][y][(z) - 1] + \
			(A)[(t) % 2][x][y][(z) + 1])

void naive_scalar(float *A, int NX, int NY, int NZ, int T);
double checksum_result(int NX, int NY, int NZ, float (*A)[NY + 2 * YSTART][NZ + 2 * ZSTART]);
