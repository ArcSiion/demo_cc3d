#include "define_512.h"

#define run_and_test(func, A)	for (x = 0; x < NX + 2 * XSTART; x++) {\
		for (y = YSTART; y < NY + YSTART; y++) {\
		A[0][x][y] = A_backup[0][x][y];\
		A[1][x][y] = A_backup[0][x][y];\
		}\
	}\
	gettimeofday(&start, 0);\
	func((int *)A, NX, NY, T);\
	gettimeofday(&end, 0);\
	if(checkresult(NX, NY, ( int (* )[ NY+2*YSTART])&A_correct[T % 2][0][0], ( int (* )[ NY+2*YSTART])&A[T % 2][0][0])) printf("Correct!\t");\
	printf(#func", NX = %d, NY = %d, T = %d, GStencil/s = %f\n", NX, NY, T,\
			((double) NX * NY * T) / (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) * 1.0e-6) / 1000000000L)


int main(int argc, char* argv[]) {

	struct timeval start, end;

	int x, y, t;
	int NX  = atoi(argv[1]);
	int NY  = atoi(argv[2]);
	int T   = atoi(argv[3]);


	int (*A)[NX + 2 * XSTART][NY+2*YSTART] = (int(*)[NX + 2 * XSTART][NY+2*YSTART])malloc(sizeof(int) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * 2);
#ifdef CHECK
	int (*A_correct)[NX + 2 * XSTART][NY+2*YSTART] = (int(*)[NX + 2 * XSTART][NY+2*YSTART])malloc(sizeof(int) * (NX + 2 * XSTART) * (NY + 2 * YSTART) * 2);
#else
	int (*A_correct)[NX + 2 * XSTART][NY+2*YSTART] = A;
#endif
	int (*A_backup)[NX + 2 * XSTART][NY+2*YSTART] = (int(*)[NX + 2 * XSTART][NY+2*YSTART])malloc(sizeof(int) * (NX + 2 * XSTART) * (NY + 2 * YSTART));



	srand(100); 
	for (x = 0; x < NX + 2 * XSTART; x++) {
		for (y = 0; y < NY + 2 * YSTART; y++) {
			A_backup[0][x][y] = INIT;
		}
	}


	run_and_test(naive_scalar, A_correct);

	run_and_test(naive_vector, A);

	run_and_test(vectime, A);

	run_and_test(vectime_extra_array, A);

	run_and_test(vectime_512, A);

	run_and_test(vectime_extra_array_512, A);


	free(A);
#ifdef CHECK
	free(A_correct);
#endif
	free(A_backup);
	return 0;
}

int checkresult( int NX, int NY, int (* A_correct)[ NY+2*YSTART], int (* A)[ NY+2*YSTART] ) {
	int correct = 1;
	int x, y;
#ifdef CHECK
	for (x= 0; x < NX + 2 * XSTART; x++) {
		for (y = 0; y < NY + 2 * YSTART; y++) {
			if(A[x][y] != A_correct[x][y]){
				printf("x = [%d], y = [%d], Correct = %d, Wrong = %d\n", x, y, A_correct[x][y], A[x][y]);
				correct = 0;
			}
		}
	}
#endif
	return correct;
}

