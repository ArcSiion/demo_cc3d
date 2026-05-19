#include "define.h"

#define reset_field(A) \
	for (x = 0; x < NX + 2 * XSTART; x++) { \
		for (y = 0; y < NY + 2 * YSTART; y++) { \
			for (z = 0; z < NZ + 2 * ZSTART; z++) { \
				(A)[0][x][y][z] = A_backup[0][x][y][z]; \
				(A)[1][x][y][z] = A_backup[0][x][y][z]; \
			} \
		} \
	}

#define run_and_report(func, A) \
	reset_field(A); \
	gettimeofday(&start, 0); \
	func((float *)(A), NX, NY, NZ, T); \
	gettimeofday(&end, 0); \
	printf(#func ", NX = %d, NY = %d, NZ = %d, T = %d, checksum = %.6e, GStencil/s = %f\n", \
		   NX, NY, NZ, T, \
		   checksum_result(NX, NY, NZ, (float (*)[NY + 2 * YSTART][NZ + 2 * ZSTART])&(A)[T % 2][0][0][0]), \
		   ((double)NX * NY * NZ * T) / \
			   (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) * 1.0e-6) / 1000000000.0)

int main(int argc, char *argv[]) {
	struct timeval start, end;
	int x, y, z;

	if (argc != 5) {
		printf("usage: %s <NX> <NY> <NZ> <T>\n", argv[0]);
		return 0;
	}

	int NX = atoi(argv[1]);
	int NY = atoi(argv[2]);
	int NZ = atoi(argv[3]);
	int T = atoi(argv[4]);

	if (NX <= 0 || NY <= 0 || NZ <= 0 || T <= 0) {
		printf("NX, NY, NZ, and T must be positive.\n");
		return 1;
	}

	float (*A)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART] =
		(float (*)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART])
			malloc(sizeof(float) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * (NZ + 2 * ZSTART) * 2);
	float (*A_backup)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART] =
		(float (*)[NX + 2 * XSTART][NY + 2 * YSTART][NZ + 2 * ZSTART])
			malloc(sizeof(float) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * (NZ + 2 * ZSTART));

	if (!A || !A_backup) {
		printf("allocation failed.\n");
		free(A);
		free(A_backup);
		return 1;
	}

	srand(100);
	for (x = 0; x < NX + 2 * XSTART; x++) {
		for (y = 0; y < NY + 2 * YSTART; y++) {
			for (z = 0; z < NZ + 2 * ZSTART; z++) {
				A_backup[0][x][y][z] = INIT;
			}
		}
	}

	run_and_report(naive_scalar, A);

	free(A);
	free(A_backup);

	return 0;
}

double checksum_result(int NX, int NY, int NZ, float (*A)[NY + 2 * YSTART][NZ + 2 * ZSTART]) {
	double checksum = 0.0;

	for (int x = XSTART; x < NX + XSTART; x++) {
		for (int y = YSTART; y < NY + YSTART; y++) {
			for (int z = ZSTART; z < NZ + ZSTART; z++) {
				checksum += (double)A[x][y][z];
			}
		}
	}

	return checksum;
}
