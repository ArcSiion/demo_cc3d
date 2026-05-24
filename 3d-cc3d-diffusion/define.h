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

#define CC3D_NUM_CELL_TYPES 256
#define CC3D_TYPE_MEDIUM 0
#define CC3D_TYPE_CELL_A 1
#define CC3D_TYPE_CELL_B 2
#define CC3D_TYPE_BLOCKED 3
#define CC3D_TRACKED_CELL_TYPES 4

typedef unsigned char cc3d_cell_type_t;

#define INIT ((float)(rand() % 1024))

#define Boundary_no_flux(A, t, NX, NY, NZ) do { \
	int _cc3d_t = (t) % 2; \
	int _cc3d_x, _cc3d_y, _cc3d_z; \
	for (_cc3d_y = YSTART; _cc3d_y < (NY) + YSTART; _cc3d_y++) { \
		for (_cc3d_z = ZSTART; _cc3d_z < (NZ) + ZSTART; _cc3d_z++) { \
			(A)[_cc3d_t][0][_cc3d_y][_cc3d_z] = (A)[_cc3d_t][XSTART][_cc3d_y][_cc3d_z]; \
			(A)[_cc3d_t][(NX) + XSTART][_cc3d_y][_cc3d_z] = (A)[_cc3d_t][(NX) + XSTART - 1][_cc3d_y][_cc3d_z]; \
		} \
	} \
	for (_cc3d_x = 0; _cc3d_x < (NX) + 2 * XSTART; _cc3d_x++) { \
		for (_cc3d_z = ZSTART; _cc3d_z < (NZ) + ZSTART; _cc3d_z++) { \
			(A)[_cc3d_t][_cc3d_x][0][_cc3d_z] = (A)[_cc3d_t][_cc3d_x][YSTART][_cc3d_z]; \
			(A)[_cc3d_t][_cc3d_x][(NY) + YSTART][_cc3d_z] = (A)[_cc3d_t][_cc3d_x][(NY) + YSTART - 1][_cc3d_z]; \
		} \
	} \
	for (_cc3d_x = 0; _cc3d_x < (NX) + 2 * XSTART; _cc3d_x++) { \
		for (_cc3d_y = 0; _cc3d_y < (NY) + 2 * YSTART; _cc3d_y++) { \
			(A)[_cc3d_t][_cc3d_x][_cc3d_y][0] = (A)[_cc3d_t][_cc3d_x][_cc3d_y][ZSTART]; \
			(A)[_cc3d_t][_cc3d_x][_cc3d_y][(NZ) + ZSTART] = (A)[_cc3d_t][_cc3d_x][_cc3d_y][(NZ) + ZSTART - 1]; \
		} \
	} \
} while (0)

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

#define Compute_scalar_cell_type(A, CELL_TYPE, DIFF_COEF, DECAY_COEF, t, x, y, z) \
	do { \
		cc3d_cell_type_t _cc3d_type = (CELL_TYPE)[x][y][z]; \
		float _cc3d_diff = (DIFF_COEF)[_cc3d_type]; \
		float _cc3d_decay = (DECAY_COEF)[_cc3d_type]; \
		float _cc3d_center = (A)[(t) % 2][x][y][z]; \
		float _cc3d_neighbor_sum = \
			(A)[(t) % 2][(x) - 1][y][z] + \
			(A)[(t) % 2][(x) + 1][y][z] + \
			(A)[(t) % 2][x][(y) - 1][z] + \
			(A)[(t) % 2][x][(y) + 1][z] + \
			(A)[(t) % 2][x][y][(z) - 1] + \
			(A)[(t) % 2][x][y][(z) + 1]; \
		(A)[((t) + 1) % 2][x][y][z] = \
			_cc3d_diff * _cc3d_neighbor_sum + \
			(1.0f - _cc3d_decay - 6.0f * _cc3d_diff) * _cc3d_center; \
	} while (0)

#define Compute_scalar_permeability(A, CELL_TYPE, DIFF_COEF, DECAY_COEF, PERM_COEF, t, x, y, z) \
    do { \
        int _cc3d_t = (t) % 2; \
        cc3d_cell_type_t _cc3d_center_type = (CELL_TYPE)[x][y][z]; \
        float _cc3d_center = (A)[_cc3d_t][x][y][z]; \
        float _cc3d_perm_center = (PERM_COEF)[_cc3d_center_type]; \
        float _cc3d_flux = \
            (_cc3d_perm_center * (PERM_COEF)[(CELL_TYPE)[(x) - 1][y][z]]) * \
                ((A)[_cc3d_t][(x) - 1][y][z] - _cc3d_center) + \
            (_cc3d_perm_center * (PERM_COEF)[(CELL_TYPE)[(x) + 1][y][z]]) * \
                ((A)[_cc3d_t][(x) + 1][y][z] - _cc3d_center) + \
            (_cc3d_perm_center * (PERM_COEF)[(CELL_TYPE)[x][(y) - 1][z]]) * \
                ((A)[_cc3d_t][x][(y) - 1][z] - _cc3d_center) + \
            (_cc3d_perm_center * (PERM_COEF)[(CELL_TYPE)[x][(y) + 1][z]]) * \
                ((A)[_cc3d_t][x][(y) + 1][z] - _cc3d_center) + \
            (_cc3d_perm_center * (PERM_COEF)[(CELL_TYPE)[x][y][(z) - 1]]) * \
                ((A)[_cc3d_t][x][y][(z) - 1] - _cc3d_center) + \
            (_cc3d_perm_center * (PERM_COEF)[(CELL_TYPE)[x][y][(z) + 1]]) * \
                ((A)[_cc3d_t][x][y][(z) + 1] - _cc3d_center); \
        (A)[((t) + 1) % 2][x][y][z] = \
            (1.0f - (DECAY_COEF)[_cc3d_center_type]) * _cc3d_center + \
            (DIFF_COEF)[_cc3d_center_type] * _cc3d_flux; \
    } while (0)

//=====================naive_scalar_edge_variable.c

#define Compute_scalar_edge_variable(A, CELL_TYPE, DIFF_COEF, DECAY_COEF, PERM_COEF, t, x, y, z) \
    do { \
        int _cc3d_t = (t) % 2; \
        cc3d_cell_type_t _type_c = (CELL_TYPE)[x][y][z]; \
        float _center = (A)[_cc3d_t][x][y][z]; \
        float _diff_c = (DIFF_COEF)[_type_c]; \
        float _decay_c = (DECAY_COEF)[_type_c]; \
        float _perm_c = (PERM_COEF)[_type_c]; \
        \
        cc3d_cell_type_t _type_l = (CELL_TYPE)[(x) - 1][y][z]; \
        cc3d_cell_type_t _type_r = (CELL_TYPE)[(x) + 1][y][z]; \
        cc3d_cell_type_t _type_d = (CELL_TYPE)[x][(y) - 1][z]; \
        cc3d_cell_type_t _type_u = (CELL_TYPE)[x][(y) + 1][z]; \
        cc3d_cell_type_t _type_b = (CELL_TYPE)[x][y][(z) - 1]; \
        cc3d_cell_type_t _type_f = (CELL_TYPE)[x][y][(z) + 1]; \
        \
        float _left  = (A)[_cc3d_t][(x) - 1][y][z]; \
        float _right = (A)[_cc3d_t][(x) + 1][y][z]; \
        float _down  = (A)[_cc3d_t][x][(y) - 1][z]; \
        float _up    = (A)[_cc3d_t][x][(y) + 1][z]; \
        float _back  = (A)[_cc3d_t][x][y][(z) - 1]; \
        float _front = (A)[_cc3d_t][x][y][(z) + 1]; \
        \
        float _w_l = _perm_c * (PERM_COEF)[_type_l] * 0.5f * (_diff_c + (DIFF_COEF)[_type_l]); \
        float _w_r = _perm_c * (PERM_COEF)[_type_r] * 0.5f * (_diff_c + (DIFF_COEF)[_type_r]); \
        float _w_d = _perm_c * (PERM_COEF)[_type_d] * 0.5f * (_diff_c + (DIFF_COEF)[_type_d]); \
        float _w_u = _perm_c * (PERM_COEF)[_type_u] * 0.5f * (_diff_c + (DIFF_COEF)[_type_u]); \
        float _w_b = _perm_c * (PERM_COEF)[_type_b] * 0.5f * (_diff_c + (DIFF_COEF)[_type_b]); \
        float _w_f = _perm_c * (PERM_COEF)[_type_f] * 0.5f * (_diff_c + (DIFF_COEF)[_type_f]); \
        \
        float _flux = \
            _w_l * (_left  - _center) + \
            _w_r * (_right - _center) + \
            _w_d * (_down  - _center) + \
            _w_u * (_up    - _center) + \
            _w_b * (_back  - _center) + \
            _w_f * (_front - _center); \
        \
        (A)[((t) + 1) % 2][x][y][z] = \
            (1.0f - _decay_c) * _center + _flux; \
    } while (0)	

void naive_scalar(float *A, int NX, int NY, int NZ, int T);
void naive_scalar_cell_type(float *A, cc3d_cell_type_t *cellType, const float *diffCoef, const float *decayCoef,
							int NX, int NY, int NZ, int T);
void naive_scalar_permeability(float *A, cc3d_cell_type_t *cellType, const float *diffCoef,
							   const float *decayCoef, const float *permeabilityCoef,
							   int NX, int NY, int NZ, int T);
double checksum_result(int NX, int NY, int NZ, float (*A)[NY + 2 * YSTART][NZ + 2 * ZSTART]);
void init_cc3d_coefficients(float *diffCoef, float *decayCoef);
void init_cc3d_permeability(float *permeabilityCoef);
void init_cell_type_field(int NX, int NY, int NZ,
						  cc3d_cell_type_t (*cellType)[NY + 2 * YSTART][NZ + 2 * ZSTART],
						  long typeCounts[CC3D_TRACKED_CELL_TYPES]);

void naive_scalar_edge_variable(float *A,
                                cc3d_cell_type_t *cellType,
                                const float *diffCoef,
                                const float *decayCoef,
                                const float *permeabilityCoef,
                                int NX, int NY, int NZ, int T);						  