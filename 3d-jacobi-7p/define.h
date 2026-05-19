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



#define Compute_one(v_center_0, v_center_1, v_center_2, a, b)  	v_center_2 = loadv_x_c_3d(y, z + (a));   \
																v_y_minus_1 = ( y == YSTART ) ? \
																					setv_3d(x, y - YSLOPE, z + (a) - 1)    \
																				:   loadv_x_c_3d(y - 1, z  + (a) - 1);  \
																v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
																					setv_3d(x, y + YSLOPE, z + (a) - 1)   \
																				:   loadv_x_c_3d(y + 1, z + (a) - 1);   \
																Compute_1vector(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
																				v_y_plus_1, \
																				loadv_x_m_3d(y, z + (a) - 1), \
																				loadv_x_p_3d(y, z + (a) - 1)); \
																Input_Output_##b(out, v_center_0, in);\
																vstore(BV0[y - YSTART][z - ZSTART + (a) - 1][0], v_center_0)

#define Compute_last_one(v_center_0, v_center_1, v_center_2, a, b)  v_center_2 = ( z > NZ + ZSTART - VECLEN - VECLEN * 2 ) ? \
																					setv_3d(x, y, z   + (a)) \
																				:   loadv_x_c_3d(y, z   + (a));   \
																	v_y_minus_1 = ( y == YSTART ) ? \
																						setv_3d(x, y - YSLOPE, z + (a) - 1)    \
																					:   loadv_x_c_3d(y - 1, z  + (a) - 1);  \
																	v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
																						setv_3d(x, y + YSLOPE, z + (a) - 1)   \
																					:   loadv_x_c_3d(y + 1, z + (a) - 1);   \
																	Compute_1vector(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
																					v_y_plus_1, \
																					loadv_x_m_3d(y, z + (a) - 1), \
																					loadv_x_p_3d(y, z + (a) - 1)); \
																	Input_Output_##b(out, v_center_0, in);\
																	vstore(BV0[y - YSTART][z - ZSTART + (a) - 1][0], v_center_0)

#if defined(simplestencil)

		#define C1 0.143

		#define SET_COFF  vec c1;\
						vallset(c1,C1)

		#define Compute_scalar(A,t,i,j,k) A[(t+1)%2][i][j][k] =  C1 * ((((((A[t%2][i][j][k] + A[t%2][i][j][k-1]) \
																	+ A[t%2][i][j][k+1]) + A[t%2][i][j-1][k]) \
																	+ A[t%2][i][j+1][k]) + A[t%2][i-1][j][k]) \
																	+ A[t%2][i+1][j][k])

		#define Compute_1vector(vcenter,vz_minus_1,vz_plus_1,vy_minus_1, vy_plus_1,vx_minus_1, vx_plus_1)	vz_minus_1 = _mm256_mul_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(\
													_mm256_add_pd(vcenter,vz_minus_1),vz_plus_1),vy_minus_1),vy_plus_1),vx_minus_1),vx_plus_1),c1);

#elif defined(heat)

		#define C0 0.125
		#define C1 -2.0
		#define SET_COFF 		__m256d vc0 = _mm256_set1_pd(C0);\
								__m256d vc1 = _mm256_set1_pd(C1)

		#define Compute_scalar(A,t,x,y,z)	A[(t+1)%2][x][y][z] =	  0.125*(A[t%2][x - 1][y][z] - 2.0 * A[t%2][x][y][z] + A[t%2][x + 1][y][z])\
																	+ 0.125*(A[t%2][x][y - 1][z] - 2.0 * A[t%2][x][y][z] + A[t%2][x][y + 1][z])\
																	+ 0.125*(A[t%2][x][y][z - 1] - 2.0 * A[t%2][x][y][z] + A[t%2][x][y][z + 1])\
																	+ A[t%2][x][y][z]

		#define Compute_1vector(vcenter,\
								vz_minus_1, vz_plus_1,\
								vy_minus_1, vy_plus_1,\
								vx_minus_1, vx_plus_1) vz_minus_1 = _mm256_add_pd( _mm256_add_pd( _mm256_add_pd(\
																_mm256_mul_pd(vc0, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vc1, vcenter), vx_minus_1), vx_plus_1)),\
																_mm256_mul_pd(vc0, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vc1, vcenter), vy_minus_1), vy_plus_1))),\
																_mm256_mul_pd(vc0, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(vc1, vcenter), vz_minus_1), vz_plus_1))),\
																vcenter)
														//vz_minus_1 = _mm256_add_pd( _mm256_add_pd( _mm256_add_pd(\
																_mm256_mul_pd(vc0, _mm256_add_pd(_mm256_fnmadd_pd(vc1, vcenter, vx_minus_1), vx_plus_1)),\
																_mm256_mul_pd(vc0, _mm256_add_pd(_mm256_fnmadd_pd(vc1, vcenter, vy_minus_1), vy_plus_1))),\
																_mm256_mul_pd(vc0, _mm256_add_pd(_mm256_fnmadd_pd(vc1, vcenter, vz_minus_1), vz_plus_1))),\
																vcenter)
#else

#define C0 0.4
#define C1 0.1

#define SET_COFF  vec vc1, vc0 ;\
				 vallset(vc1,C1);\
				 vallset(vc0,C0)

#define Compute_scalar(A,t,i,j,k) A[(t+1)%2][i][j][k] =  C0 * A[t%2][i][j][k] + C1 * (((((( (A[t%2][i][j][k-1]) \
															 + A[t%2][i][j][k+1])) + A[t%2][i][j-1][k]) \
															 + A[t%2][i][j+1][k]) + A[t%2][i-1][j][k]) \
															 + A[t%2][i+1][j][k])

#define Compute_1vector(vcenter,vz_minus_1,vz_plus_1,vy_minus_1, vy_plus_1,vx_minus_1, vx_plus_1)	vz_minus_1 = _mm256_fmadd_pd(vcenter, vc0,\
                                                                                        _mm256_mul_pd(\
                                                                                        _mm256_add_pd(\
                                                                                        _mm256_add_pd(\
                                                                                        _mm256_add_pd(\
                                                                                        _mm256_add_pd(\
                                                                                        _mm256_add_pd(vz_minus_1,vz_plus_1),vy_minus_1),vy_plus_1),vx_minus_1),vx_plus_1),vc1));

#endif


void naive_scalar(double* A, int NX, int NY, int NZ, int T);
void naive_vector(double* A, int NX, int NY, int NZ, int T);
void vectime(double* A, int NX, int NY, int NZ, int T);
void vectime_extra_array(double* A, int NX, int NY, int NZ, int T);
void vectime_extra_array_unroll8(double* A, int NX, int NY, int NZ, int T);
void vectime_transpose_boundary(double* A, int NX, int NY, int NZ, int T);
void vectime_transpose_boundary_extra_array(double* A, int NX, int NY, int NZ, int T);
int checkresult( int NX, int NY, int NZ, double (* A_correct)[ NY+2*YSTART][NZ+2*ZSTART], double (* A)[ NY+2*YSTART][NZ+2*ZSTART]);
