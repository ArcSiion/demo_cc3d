
#include "define.h"

#ifndef __AVX512F__
#error "define_512.h requires AVX-512F. Compile with -mavx512f."
#endif

typedef __m512d vec_512;

#define back_scalar_512(NX, NY, NZ, T) \
    ((T) >= VVECLEN && (NX) >= VVECLEN && (NY) >= VVECLEN && (NZ) >= VVECLEN && \
     (NX) < STRIDE * (VVECLEN - 1) + 1)

#define vload_512(a,b) (a) = _mm512_loadu_pd(&(b))
#define vstore_512(a,b) _mm512_storeu_pd(&(a),(b))
#define vloadset_512(a, h, g, f, e, d, c, b, in_a) { a = _mm512_set_pd(h, g, f, e, d, c, b, in_a); }
#define vallset_512(a,b) { a = _mm512_set1_pd(b); }

#define alloc_extra_array_512(size) ({  \
        void *ptr = NULL;                       \
        if (posix_memalign(&ptr, 64, size) != 0) { ptr = NULL; } \
        ptr;                                \
})

#define free_extra_array_512(a) free_extra_array(a)

#define load_v2_512(x, y, z, xshift, yshift, zshift) setv_3d_512((x) + (xshift), (y) + (yshift), (z) + (zshift))

#define Add_4_vectors_512(a, b, c, d) _mm512_add_pd(\
                                    _mm512_add_pd(\
                                    _mm512_add_pd(a, b), c), d)

#ifdef simplestencil

#define SET_COFF_512  vec_512 vc1_512;\
                      vallset_512(vc1_512,C1)

#define Compute_1vector_512(v_center_0, \
                            v_center_1, \
                            v_center_2, \
                            v_all_d_1_0, \
                            v_all_d_1_1, \
                            v_all_d_1_2, \
                            v_all_d_2_0, \
                            v_all_d_2_1, \
                            v_all_d_2_2) v_center_0 = _mm512_mul_pd(\
                                                            _mm512_add_pd(\
                                                            _mm512_add_pd(\
                                                            _mm512_add_pd(\
                                                            _mm512_add_pd(\
                                                            _mm512_add_pd(\
                                                            _mm512_add_pd(\
                                                            _mm512_add_pd(\
                                                            _mm512_add_pd(v_center_0, \
                                                                          v_center_1), \
                                                                          v_center_2), \
                                                                          v_all_d_1_0), \
                                                                          v_all_d_1_1), \
                                                                          v_all_d_1_2), \
                                                                          v_all_d_2_0), \
                                                                          v_all_d_2_1), \
                                                                          v_all_d_2_2), \
                                                                          vc1_512)

#else

#define SET_COFF_512  vec_512 vc1_512, vc0_512, vc2_512, vc3_512 ;\
                      vallset_512(vc1_512,C1);\
                      vallset_512(vc0_512,C0);\
                      vallset_512(vc2_512,C2);\
                      vallset_512(vc3_512,C3)

#define Compute_1vector_512(v_center_0, \
                            v_center_1, \
                            v_center_2, \
                            v_all_d_1_0, \
                            v_all_d_1_1, \
                            v_all_d_1_2, \
                            v_all_d_2_0, \
                            v_all_d_2_1, \
                            v_all_d_2_2) v_center_0 = _mm512_fmadd_pd(vc3_512,\
                                                            _mm512_add_pd(_mm512_add_pd(v_all_d_2_1,v_all_d_1_0),v_all_d_1_2), \
                                                            _mm512_fmadd_pd(vc1_512,\
                                                                            _mm512_add_pd(_mm512_add_pd(v_center_0,v_all_d_1_1),v_center_2), \
                                                                            _mm512_fmadd_pd(vc2_512,\
                                                                                            _mm512_add_pd(v_all_d_2_0, v_all_d_2_2),\
                                                                                            _mm512_mul_pd(vc0_512,v_center_1))))

#endif

#define vrotate_512_low2high(a) \
    _mm512_castsi512_pd(_mm512_alignr_epi64(_mm512_castpd_si512(a), _mm512_castpd_si512(a), 1))

#define vrotate_512_high2low(a) \
    _mm512_castsi512_pd(_mm512_alignr_epi64(_mm512_castpd_si512(a), _mm512_castpd_si512(a), 7))

#define transpose8x8_pd_512(r0,r1,r2,r3,r4,r5,r6,r7) do {                       \
    __m512d t0  = _mm512_unpacklo_pd((r0),(r1));                                \
    __m512d t1  = _mm512_unpackhi_pd((r0),(r1));                                \
    __m512d t2  = _mm512_unpacklo_pd((r2),(r3));                                \
    __m512d t3  = _mm512_unpackhi_pd((r2),(r3));                                \
    __m512d t4  = _mm512_unpacklo_pd((r4),(r5));                                \
    __m512d t5  = _mm512_unpackhi_pd((r4),(r5));                                \
    __m512d t6  = _mm512_unpacklo_pd((r6),(r7));                                \
    __m512d t7  = _mm512_unpackhi_pd((r6),(r7));                                \
    const __m512i idx0 = _mm512_setr_epi64(0, 1, 8, 9, 4, 5, 12, 13);           \
    const __m512i idx1 = _mm512_setr_epi64(2, 3,10,11, 6, 7, 14, 15);           \
    __m512d tt0 = _mm512_permutex2var_pd(t0, idx0, t2);                         \
    __m512d tt1 = _mm512_permutex2var_pd(t0, idx1, t2);                         \
    __m512d tt2 = _mm512_permutex2var_pd(t1, idx0, t3);                         \
    __m512d tt3 = _mm512_permutex2var_pd(t1, idx1, t3);                         \
    __m512d tt4 = _mm512_permutex2var_pd(t4, idx0, t6);                         \
    __m512d tt5 = _mm512_permutex2var_pd(t4, idx1, t6);                         \
    __m512d tt6 = _mm512_permutex2var_pd(t5, idx0, t7);                         \
    __m512d tt7 = _mm512_permutex2var_pd(t5, idx1, t7);                         \
    (r0) = _mm512_shuffle_f64x2(tt0, tt4, 0x44);                                \
    (r1) = _mm512_shuffle_f64x2(tt2, tt6, 0x44);                                \
    (r2) = _mm512_shuffle_f64x2(tt1, tt5, 0x44);                                \
    (r3) = _mm512_shuffle_f64x2(tt3, tt7, 0x44);                                \
    (r4) = _mm512_shuffle_f64x2(tt0, tt4, 0xEE);                                \
    (r5) = _mm512_shuffle_f64x2(tt2, tt6, 0xEE);                                \
    (r6) = _mm512_shuffle_f64x2(tt1, tt5, 0xEE);                                \
    (r7) = _mm512_shuffle_f64x2(tt3, tt7, 0xEE);                                \
} while(0)

#define transpose_512(a,b,c,d,e,f,g,h, in, out) do { \
    transpose8x8_pd_512(a,b,c,d,e,f,g,h); \
} while(0)

#define Input_Output_1_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = (v1); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
    (in) = vrotate_512_low2high(in); \
} while(0)

#define Input_Output_2_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = _mm512_mask_blend_pd(0x02, (out), _mm512_broadcastsd_pd(_mm512_castpd512_pd128(v1))); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
    (in) = vrotate_512_low2high(in); \
} while(0)

#define Input_Output_3_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = _mm512_mask_blend_pd(0x04, (out), _mm512_broadcastsd_pd(_mm512_castpd512_pd128(v1))); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
    (in) = vrotate_512_low2high(in); \
} while(0)

#define Input_Output_4_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = _mm512_mask_blend_pd(0x08, (out), _mm512_broadcastsd_pd(_mm512_castpd512_pd128(v1))); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
    (in) = vrotate_512_low2high(in); \
} while(0)

#define Input_Output_5_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = _mm512_mask_blend_pd(0x10, (out), _mm512_broadcastsd_pd(_mm512_castpd512_pd128(v1))); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
    (in) = vrotate_512_low2high(in); \
} while(0)

#define Input_Output_6_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = _mm512_mask_blend_pd(0x20, (out), _mm512_broadcastsd_pd(_mm512_castpd512_pd128(v1))); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
    (in) = vrotate_512_low2high(in); \
} while(0)

#define Input_Output_7_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = _mm512_mask_blend_pd(0x40, (out), _mm512_broadcastsd_pd(_mm512_castpd512_pd128(v1))); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
    (in) = vrotate_512_low2high(in); \
} while(0)

#define Input_Output_8_512(out, v1, in) do { \
    (v1) = vrotate_512_high2low(v1); \
    (out) = _mm512_mask_blend_pd(0x80, (out), _mm512_broadcastsd_pd(_mm512_castpd512_pd128(v1))); \
    (v1) = _mm512_mask_blend_pd(0x01, (v1), (in)); \
} while(0)

#define setv_3d_512(x, y, z) \
    _mm512_set_pd(B[(t+1)%2][(x)           ][(y)][(z)], \
                  B[(t)%2  ][(x) + STRIDE  ][(y)][(z)], \
                  B[(t+1)%2][(x) + STRIDE*2][(y)][(z)], \
                  B[(t)%2  ][(x) + STRIDE*3][(y)][(z)], \
                  B[(t+1)%2][(x) + STRIDE*4][(y)][(z)], \
                  B[(t)%2  ][(x) + STRIDE*5][(y)][(z)], \
                  B[(t+1)%2][(x) + STRIDE*6][(y)][(z)], \
                  B[(t)%2  ][(x) + STRIDE*7][(y)][(z)])

#define vset_3d_512(vec, x, y, z) (vec) = setv_3d_512((x), (y), (z))

#define load_x_m_512(x, y, z, yshift, zshift) _mm512_loadu_pd(&(BV0[(y) + (yshift)][(z) + (zshift) - ZSTART][0]))
#define load_x_c_512(x, y, z, yshift, zshift) _mm512_loadu_pd(&(BV1[(y) + (yshift)][(z) + (zshift) - ZSTART][0]))
#define load_x_p_512(x, y, z, yshift, zshift) _mm512_loadu_pd(&(BV2[(y) + (yshift)][(z) + (zshift) - ZSTART][0]))

void vectime_transpose_boundary_extra_array_512(double* A, int NX, int NY, int NZ, int T);

#define vstore_set_3d_512(x, y, z, v) do { \
    _mm512_storeu_pd(tmp, (v)); \
    B[(t + 1) % 2][(x) + STRIDE * 7][(y)][(z)] = tmp[0]; \
    B[(t)     % 2][(x) + STRIDE * 6][(y)][(z)] = tmp[1]; \
    B[(t + 1) % 2][(x) + STRIDE * 5][(y)][(z)] = tmp[2]; \
    B[(t)     % 2][(x) + STRIDE * 4][(y)][(z)] = tmp[3]; \
    B[(t + 1) % 2][(x) + STRIDE * 3][(y)][(z)] = tmp[4]; \
    B[(t)     % 2][(x) + STRIDE * 2][(y)][(z)] = tmp[5]; \
    B[(t + 1) % 2][(x) + STRIDE * 1][(y)][(z)] = tmp[6]; \
    B[(t)     % 2][(x) + STRIDE * 0][(y)][(z)] = tmp[7]; \
} while (0)
