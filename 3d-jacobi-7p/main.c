#include "define_512.h"

#define run_and_test(func, A)	for (x = 0; x < NX + 2 * XSTART; x++) {\
		for (y = 0; y < NY + 2 * YSTART; y++) {\
			for (z = 0; z < NZ + 2 * ZSTART; z++) {\
				A[0][x][y][z] = A_backup[0][x][y][z];\
				A[1][x][y][z] = A_backup[0][x][y][z];\
			}\
		}\
	}\
	gettimeofday(&start, 0);\
	func((double *)A, NX, NY, NZ, T);\
	gettimeofday(&end, 0);\
	if(checkresult(NX, NY, NZ, ( double (* )[NY+2*YSTART][ NZ+2*ZSTART])&A_correct[T % 2][0][0][0], ( double (* )[NY+2*YSTART][ NZ+2*ZSTART])&A[T % 2][0][0][0])) printf("Correct!\t");\
	printf(#func", NX = %d, NY = %d, NZ = %d, T = %d, GStencil/s = %f\n", \
	 								NX, NY, NZ, T, ((double) NX * NY * NZ * T) / (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) * 1.0e-6) / 1000000000L)


int main(int argc, char* argv[]) {

	struct timeval start, end;
	long int i, j, k;
	int x, y, z, t;
	if (argc != 5) {
		printf("usage: %s <NX> <NY> <NZ> <T>\n", argv[0]);
		return 0;
	}
	int NX = atoi(argv[1]);
	int NY = atoi(argv[2]);
	int NZ = atoi(argv[3]);
	int T = atoi(argv[4]);
	
	double (*A)[NX + 2 * XSTART][NY+2*YSTART][NZ+2*ZSTART] = (double(*)[NX + 2 * XSTART][NY+2*YSTART][NZ+2*ZSTART])malloc(sizeof(double) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * (NZ + 2 * ZSTART) * 2);
#ifdef CHECK
	double (*A_correct)[NX + 2 * XSTART][NY+2*YSTART][NZ+2*ZSTART] = (double(*)[NX + 2 * XSTART][NY+2*YSTART][NZ+2*ZSTART])malloc(sizeof(double) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * (NZ + 2 * ZSTART) * 2);
#else
	double (*A_correct)[NX + 2 * XSTART][NY+2*YSTART][NZ+2*ZSTART] = A;
#endif
	double (*A_backup)[NX + 2 * XSTART][NY+2*YSTART][NZ+2*ZSTART] = (double(*)[NX + 2 * XSTART][NY+2*YSTART][NZ+2*ZSTART])malloc(sizeof(double) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * (NZ + 2 * ZSTART) ) ;

	srand(100);
	for (i = 0; i < NX + 2 * XSTART; i++) {
		for (j = 0; j < NY + 2 * YSTART; j++) {
			for (k = 0; k < NZ + 2 * ZSTART; k++) {
				A_backup[0][i][j][k] = INIT;
			}
		}
	}

	run_and_test(naive_scalar, A_correct);
	run_and_test(naive_vector, A);
	run_and_test(vectime, A);
	run_and_test(vectime_extra_array, A);
	run_and_test(vectime_extra_array_unroll8, A);
	run_and_test(vectime_transpose_boundary, A);
	run_and_test(vectime_transpose_boundary_extra_array, A);
	run_and_test(vectime_512, A);
	run_and_test(vectime_512_hybrid, A);
	//run_and_test(vectime_512_debug, A);
	run_and_test(vectime_extra_array_512, A);
	run_and_test(vectime_extra_array_512_hybrid, A);
	run_and_test(vectime_extra_array_unroll8_512, A);
	run_and_test(vectime_extra_array_unroll8_512_hybrid, A);
	run_and_test(vectime_transpose_boundary_512, A);
	run_and_test(vectime_transpose_boundary_512_hybrid, A);
	run_and_test(vectime_transpose_boundary_extra_array_512, A);
	run_and_test(vectime_transpose_boundary_extra_array_512_hybrid, A);
	
#ifdef CHECK
	free(A_correct);
#endif
	free(A);
	free(A_backup);	

	return 0;
}

int checkresult( int NX, int NY, int NZ, double (* A_correct)[NY+2*YSTART][ NZ+2*ZSTART], double (* A)[NY+2*YSTART][ NZ+2*ZSTART] ) {
	int correct = 1;
#ifdef CHECK
	int x, y, z;
	for (x= XSTART - XSLOPE; x < NX + XSTART + XSLOPE; x++) {
		for (y = YSTART - YSLOPE; y < NY + YSTART + YSLOPE; y++) {
			for (z = ZSTART - ZSLOPE; z < NZ + ZSTART + ZSLOPE; z++) {
	// for (x= XSTART; x < NX + XSTART; x++) {
	// 	for (y = YSTART; y < NY + YSTART; y++) {
	// 		for (z = ZSTART; z < NZ + ZSTART; z++) {
				if(A[x][y][z] != A_correct[x][y][z]){
					printf("x = [%d], y = [%d], z = [%d], Correct = %f, Wrong = %f\n", x-XSTART+1, y-YSTART+1, z, A_correct[x][y][z], A[x][y][z]);
					correct = 0;
				}
			}
		}
	}
#endif
	return correct;
}