#include "../common.h"

#define STRIDE 2

#define XSTART 1
#define YSTART 4
#define YSLOPE 1
#define XSLOPE 1

#define INIT (rand() % 2)

#define VECLEN_INT 8


#define Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                        v_x_minus_1, v_center_1, v_x_plus_1,\
                        v_x_minus_0, v_center_0, v_x_plus_0)  v_center_0 =  _mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(\
						v_x_minus_0, v_center_0),\
                        v_x_minus_1), v_x_plus_0), v_x_plus_1),\
                        v_x_minus_2), v_center_2), v_x_plus_2);\
						v_x_plus_0 = _mm256_blendv_epi8(vone, vzero, _mm256_or_si256(_mm256_cmpgt_epi32(vtwo, v_center_0),_mm256_cmpgt_epi32( v_center_0, vthree)));\
						v_x_minus_0 = _mm256_blendv_epi8(v_center_1, vone, _mm256_cmpeq_epi32(v_center_0, vthree));\
						v_center_0 = _mm256_blendv_epi8(v_x_minus_0, v_x_plus_0, _mm256_cmpeq_epi32(v_center_1, vone))


						


#define  Compute_scalar(A, t, x, y)  A[(t+1)%2][x][y] = 	A[t%2][x-1][y+1] + A[t%2][x-1][y] + \
															A[t%2][x-1][y-1] + A[t%2][x][y+1] + \
															A[t%2][x][y-1]   + A[t%2][x+1][y+1] + \
															A[t%2][x+1][y]   + A[t%2][x+1][y-1];\
															if(A[(t)%2][x][y] == 1) {\
																if((A[(t+1)%2][x][y] < 2) || (A[(t+1)%2][x][y] > 3)){\
																	A[(t+1)%2][x][y] = 0;\
																} else {\
																	A[(t+1)%2][x][y] = 1;\
																}\
															} else { \
																if( A[(t+1)%2][x][y] == 3) {\
																	A[(t+1)%2][x][y] = 1;\
																} else{ \
																	A[(t+1)%2][x][y] = A[(t)%2][x][y];\
																}\
															}
															


void vectime(int* A, int NX, int NY, int T);
void vectime_extra_array(int* A, int NX, int NY, int T);
void naive_scalar(int * A, int NX, int NY, int T);
void naive_vector(int * A, int NX, int NY, int T);
int checkresult( int NX, int NY, int (* A_correct)[ NY+2*YSTART], int (* A)[ NY+2*YSTART]);


