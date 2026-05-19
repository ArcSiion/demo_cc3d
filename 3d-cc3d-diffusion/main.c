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

#define run_and_report_cell_type(func, A, CELL_TYPE, DIFF_COEF, DECAY_COEF) \
	reset_field(A); \
	gettimeofday(&start, 0); \
	func((float *)(A), (cc3d_cell_type_t *)(CELL_TYPE), DIFF_COEF, DECAY_COEF, NX, NY, NZ, T); \
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
	cc3d_cell_type_t (*cellType)[NY + 2 * YSTART][NZ + 2 * ZSTART] =
		(cc3d_cell_type_t (*)[NY + 2 * YSTART][NZ + 2 * ZSTART])
			malloc(sizeof(cc3d_cell_type_t) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * (NZ + 2 * ZSTART));
	float diffCoef[CC3D_NUM_CELL_TYPES];
	float decayCoef[CC3D_NUM_CELL_TYPES];
	long typeCounts[3] = {0, 0, 0};

	if (!A || !A_backup || !cellType) {
		printf("allocation failed.\n");
		free(A);
		free(A_backup);
		free(cellType);
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

	init_cc3d_coefficients(diffCoef, decayCoef);
	init_cell_type_field(NX, NY, NZ, cellType, typeCounts);
	printf("cell_type_counts, Medium = %ld, CellA = %ld, CellB = %ld\n",
		   typeCounts[CC3D_TYPE_MEDIUM], typeCounts[CC3D_TYPE_CELL_A], typeCounts[CC3D_TYPE_CELL_B]);

	run_and_report(naive_scalar, A);
	run_and_report_cell_type(naive_scalar_cell_type, A, cellType, diffCoef, decayCoef);

	free(A);
	free(A_backup);
	free(cellType);

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

void init_cc3d_coefficients(float *diffCoef, float *decayCoef) {
	for (int i = 0; i < CC3D_NUM_CELL_TYPES; i++) {
		diffCoef[i] = CC3D_DIFFUSION_COEF;
		decayCoef[i] = CC3D_DECAY_COEF;
	}

	diffCoef[CC3D_TYPE_MEDIUM] = 0.10f;
	decayCoef[CC3D_TYPE_MEDIUM] = 0.00f;
	diffCoef[CC3D_TYPE_CELL_A] = 0.06f;
	decayCoef[CC3D_TYPE_CELL_A] = 0.02f;
	diffCoef[CC3D_TYPE_CELL_B] = 0.12f;
	decayCoef[CC3D_TYPE_CELL_B] = 0.01f;
}

void init_cell_type_field(int NX, int NY, int NZ,
						  cc3d_cell_type_t (*cellType)[NY + 2 * YSTART][NZ + 2 * ZSTART],
						  long typeCounts[3]) {
	for (int x = 0; x < NX + 2 * XSTART; x++) {
		for (int y = 0; y < NY + 2 * YSTART; y++) {
			for (int z = 0; z < NZ + 2 * ZSTART; z++) {
				cc3d_cell_type_t type = CC3D_TYPE_MEDIUM;

				if (x >= XSTART && x < NX + XSTART &&
					y >= YSTART && y < NY + YSTART &&
					z >= ZSTART && z < NZ + ZSTART) {
					int relX = x - XSTART;
					int relY = y - YSTART;
					int relZ = z - ZSTART;

					if (relX >= NX / 4 && relX < (3 * NX) / 4 &&
						relY >= NY / 4 && relY < (3 * NY) / 4 &&
						relZ >= NZ / 4 && relZ < (3 * NZ) / 4) {
						type = CC3D_TYPE_CELL_A;
					} else if (((relX / 4) + (relY / 4) + (relZ / 4)) % 5 == 0) {
						type = CC3D_TYPE_CELL_B;
					}

					typeCounts[type]++;
				}

				cellType[x][y][z] = type;
			}
		}
	}
}
