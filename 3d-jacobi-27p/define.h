#include "../common.h"

typedef __m256d vec;

#define XSTART 1
#define YSTART 1
#define ZSTART 4

#define XSLOPE 1
#define YSLOPE 1
#define ZSLOPE 1

#define STRIDE 2

#define INIT 1.0 * (rand() % 1024)

#define back_scalar_256(NX, NY, NZ, T) \
	((T) >= VECLEN && (NX) >= VECLEN && (NY) >= VECLEN && (NZ) >= VECLEN && \
	 (NX) < STRIDE * (VECLEN - 1) + 1)




#define	load_v2(x, y, z, xshift, yshift, zshift)	setv_3d((x) + (xshift), (y) + (yshift), (z) + (zshift))


#define Add_4_vectors(a, b, c, d)	_mm256_add_pd(\
									_mm256_add_pd(\
									_mm256_add_pd(a, b), c), d)\



#ifdef simplestencil


#define C1 0.037


#define SET_COFF  __m256d vc1;\
				 vallset(vc1,C1);

#define Compute_scalar(A,t,i,j,k) A[(t+1)%2][i][j][k] =  C1 * (((((((((\
																A[t%2][i]	[j]		[k]  \
															+	A[t%2][i-1]	[j]		[k]) \
															+ 	A[t%2][i+1]	[j]		[k]) \
															+	A[t%2][i]	[j-1]	[k]) \
															+ 	A[t%2][i-1]	[j-1]	[k]) \
															+	A[t%2][i+1]	[j-1]	[k]) \
															+ 	A[t%2][i]	[j+1]	[k]) \
															+	A[t%2][i-1]	[j+1]	[k]) \
															+ 	A[t%2][i+1]	[j+1]	[k]) \
															+	((((((((A[t%2][i]	[j]		[k-1] \
															+	A[t%2][i-1]	[j]		[k-1]) \
															+	A[t%2][i+1]	[j]		[k-1]) \
															+	A[t%2][i]	[j-1]	[k-1]) \
															+	A[t%2][i-1]	[j-1]	[k-1]) \
															+	A[t%2][i+1]	[j-1]	[k-1]) \
															+	A[t%2][i]	[j+1]	[k-1]) \
															+	A[t%2][i-1]	[j+1]	[k-1]) \
															+	A[t%2][i+1]	[j+1]	[k-1]) \
															+ 	((((((((A[t%2][i]	[j]		[k+1] \
															+	A[t%2][i-1]	[j]		[k+1]) \
															+	A[t%2][i+1]	[j]		[k+1]) \
															+	A[t%2][i]	[j-1]	[k+1]) \
															+	A[t%2][i-1]	[j-1]	[k+1]) \
															+	A[t%2][i+1]	[j-1]	[k+1]) \
															+	A[t%2][i]	[j+1]	[k+1]) \
															+	A[t%2][i-1]	[j+1]	[k+1]) \
															+	A[t%2][i+1]	[j+1]	[k+1]))


#define 			Compute_1vector(v_center_0, \
                                    v_center_1, \
                                    v_center_2, \
                                    v_all_d_1_0, \
                                    v_all_d_1_1, \
                                    v_all_d_1_2, \
                                    v_all_d_2_0, \
                                    v_all_d_2_1, \
                                    v_all_d_2_2) 	v_center_0 = _mm256_mul_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(	v_center_0, \
																	v_center_1), \
																	v_center_2), \
																	v_all_d_1_0), \
																	v_all_d_1_1), \
																	v_all_d_1_2), \
																	v_all_d_2_0), \
																	v_all_d_2_1), \
																	v_all_d_2_2), \
																	vc1)


#else

#define C0 0.54
#define C1 0.03
#define C2 0.02
#define C3 0.01

#define SET_COFF  vec vc1, vc0, vc2, vc3 ;\
				 vallset(vc1,C1);\
				 vallset(vc0,C0);\
				 vallset(vc2,C2);\
				 vallset(vc3,C3)

#define Compute_scalar(A,t,i,j,k)  A[(t+1)%2][i][j][k] = C0 * A[(t)%2][i][j][k] + \
                                        C1 * (A[(t)%2][i-1][j][k] + A[(t)%2][i+1][j][k] + A[(t)%2][i][j-1][k] +\
                                                A[(t)%2][i][j+1][k] + A[(t)%2][i][j][k-1] + A[(t)%2][i][j][k+1]) +\
                                        C2 * (A[(t)%2][i-1][j-1][k-1] + A[(t)%2][i-1][j-1][k+1] + A[(t)%2][i-1][j+1][k-1] +\
                                                A[(t)%2][i-1][j+1][k+1] + A[(t)%2][i+1][j-1][k-1] + A[(t)%2][i+1][j-1][k+1] +\
                                                A[(t)%2][i+1][j+1][k-1] + A[(t)%2][i+1][j+1][k+1]) + \
                                        C3 * (A[(t)%2][i-1][j-1][k] + A[(t)%2][i-1][j][k-1] + A[(t)%2][i-1][j][k+1] +\
                                                A[(t)%2][i-1][j+1][k] + A[(t)%2][i][j-1][k-1] + A[(t)%2][i][j-1][k+1] +\
                                                A[(t)%2][i][j+1][k-1] + A[(t)%2][i][j+1][k+1] + A[(t)%2][i+1][j-1][k] +\
                                                A[(t)%2][i+1][j][k-1] + A[(t)%2][i+1][j][k+1] + A[(t)%2][i+1][j+1][k])


#define 			Compute_1vector(v_center_0, \
                                    v_center_1, \
                                    v_center_2, \
                                    v_all_d_1_0, \
                                    v_all_d_1_1, \
                                    v_all_d_1_2, \
                                    v_all_d_2_0, \
                                    v_all_d_2_1, \
                                    v_all_d_2_2) 	v_center_0 = _mm256_fmadd_pd(vc3,\
																	_mm256_add_pd(_mm256_add_pd(v_all_d_2_1,v_all_d_1_0),v_all_d_1_2), \
																	_mm256_fmadd_pd(vc1,\
																					_mm256_add_pd(_mm256_add_pd(v_center_0,v_all_d_1_1),v_center_2), \
																					_mm256_fmadd_pd(vc2,\
																									_mm256_add_pd(v_all_d_2_0, v_all_d_2_2),\
																									_mm256_mul_pd(vc0,v_center_1))))

#endif



void naive_scalar(double* A, int NX, int NY, int NZ, int T);
void naive_vector(double* A, int NX, int NY, int NZ, int T);
void vectime(double* A, int NX, int NY, int NZ, int T);
void vectime_extra_array(double* A, int NX, int NY, int NZ, int T);
void vectime_transpose_boundary(double* A, int NX, int NY, int NZ, int T);
void vectime_transpose_boundary_extra_array(double* A, int NX, int NY, int NZ, int T);
void vectime_transpose_boundary_extra_array_512(double* A, int NX, int NY, int NZ, int T);
void vectime3d_2(double* A, int NX, int NY, int NZ, int T);
int checkresult( int NX, int NY, int NZ, double (* A_correct)[ NY+2*YSTART][NZ+2*ZSTART], double (* A)[ NY+2*YSTART][NZ+2*ZSTART]);
