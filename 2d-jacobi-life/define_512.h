#include "define.h"

#define VECLEN_INT_512 16

typedef __m512i veci_512;


#ifdef __linux__
#define alloc_extra_array_512(a) _mm_malloc(a, 64)
#define free_extra_array_512(a) _mm_free(a)
#else
#define alloc_extra_array_512(a) malloc(a)
#define free_extra_array_512(a) free(a)
#endif


#define vloadi_512(a) _mm512_loadu_si512((__m512i * )&a)
#define vloadi2_512(b,a) b = _mm512_loadu_si512((__m512i * )&a)
#define vstorei_512(a,b) _mm512_storeu_si512((__m512i * ) &a, b)

#define vallseti_512(a,b) a = _mm512_set1_epi32(b)

#define SET_LIFE_CONST_512      vzero = _mm512_set1_epi32(0);\
                                vone = _mm512_set1_epi32(1);\
                                vtwo = _mm512_set1_epi32(2);\
                                vthree = _mm512_set1_epi32(3)


#define vrotate_high2lowi_512(a)      _mm512_alignr_epi32(a, a, 15)
#define vrotate_low2highi_512(a)      _mm512_alignr_epi32(a, a, 1)

#define vloadseti_512(v, B, t, x, y)     v = _mm512_set_epi32( B[ (t + 1)%2 ][ x                 ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 1    ][ y ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 2    ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 3    ][ y ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 4    ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 5    ][ y ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 6    ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 7    ][ y ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 8    ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 9    ][ y ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 10   ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 11   ][ y ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 12   ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 13   ][ y ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 14   ][ y ], \
                                                               B[ (t)%2     ][ x + STRIDE * 15   ][ y ])

#define vloadseti_blk_512(v, B, t, x, y) v = _mm512_set_epi32( B[ (t + 1)%2 ][ x                 ][ y + 15 ], \
                                                               B[ (t)%2     ][ x + STRIDE * 1    ][ y + 14 ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 2    ][ y + 13 ], \
                                                               B[ (t)%2     ][ x + STRIDE * 3    ][ y + 12 ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 4    ][ y + 11 ], \
                                                               B[ (t)%2     ][ x + STRIDE * 5    ][ y + 10 ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 6    ][ y + 9  ], \
                                                               B[ (t)%2     ][ x + STRIDE * 7    ][ y + 8  ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 8    ][ y + 7  ], \
                                                               B[ (t)%2     ][ x + STRIDE * 9    ][ y + 6  ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 10   ][ y + 5  ], \
                                                               B[ (t)%2     ][ x + STRIDE * 11   ][ y + 4  ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 12   ][ y + 3  ], \
                                                               B[ (t)%2     ][ x + STRIDE * 13   ][ y + 2  ], \
                                                               B[ (t + 1)%2 ][ x + STRIDE * 14   ][ y + 1  ], \
                                                               B[ (t)%2     ][ x + STRIDE * 15   ][ y      ])


#define Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,\
                            v_x_minus_1, v_center_1, v_x_plus_1,\
                            v_x_minus_0, v_center_0, v_x_plus_0)  \
                            do {\
                                v_center_0 = _mm512_add_epi32(_mm512_add_epi32(_mm512_add_epi32(_mm512_add_epi32(_mm512_add_epi32(_mm512_add_epi32(_mm512_add_epi32(\
                                             v_x_minus_0, v_center_0),\
                                             v_x_minus_1), v_x_plus_0), v_x_plus_1),\
                                             v_x_minus_2), v_center_2), v_x_plus_2);\
                                __mmask16 life_birth_mask_512 = _mm512_cmp_epi32_mask(v_center_0, vthree, _MM_CMPINT_EQ);\
                                __mmask16 life_survive_mask_512 = _mm512_cmp_epi32_mask(v_center_0, vtwo, _MM_CMPINT_EQ) & _mm512_cmp_epi32_mask(v_center_1, vone, _MM_CMPINT_EQ);\
                                v_center_0 = _mm512_maskz_mov_epi32(life_birth_mask_512 | life_survive_mask_512, vone);\
                            } while(0)



// Each Output_i appends lane 0 of v1 to the high end of out; after 16 appends,
// the previous contents have shifted out and lanes 0..15 hold the final vector.
#define Output_i_1_512(out, v1)      out = _mm512_alignr_epi32(v1, _mm512_setzero_si512(), 1)
#define Output_i_2_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_3_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_4_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_5_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_6_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_7_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_8_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_9_512(out, v1)      out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_10_512(out, v1)     out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_11_512(out, v1)     out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_12_512(out, v1)     out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_13_512(out, v1)     out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_14_512(out, v1)     out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_15_512(out, v1)     out = _mm512_alignr_epi32(v1, out, 1)
#define Output_i_16_512(out, v1)     out = _mm512_alignr_epi32(v1, out, 1)


#define Input_i_1_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_2_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_3_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_4_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_5_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_6_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_7_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_8_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_9_512(v1, in)        v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_10_512(v1, in)       v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_11_512(v1, in)       v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_12_512(v1, in)       v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_13_512(v1, in)       v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_14_512(v1, in)       v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_15_512(v1, in)       v1 = _mm512_mask_blend_epi32(0x0001, v1, in);\
                                     in = vrotate_low2highi_512(in)
#define Input_i_16_512(v1, in)       v1 = _mm512_mask_blend_epi32(0x0001, v1, in)


#define Input_Output_i_1_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_1_512(out, v1);\
                                             Input_i_1_512(v1, in)

#define Input_Output_i_2_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_2_512(out, v1);\
                                             Input_i_2_512(v1, in)

#define Input_Output_i_3_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_3_512(out, v1);\
                                             Input_i_3_512(v1, in)

#define Input_Output_i_4_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_4_512(out, v1);\
                                             Input_i_4_512(v1, in)

#define Input_Output_i_5_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_5_512(out, v1);\
                                             Input_i_5_512(v1, in)

#define Input_Output_i_6_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_6_512(out, v1);\
                                             Input_i_6_512(v1, in)

#define Input_Output_i_7_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_7_512(out, v1);\
                                             Input_i_7_512(v1, in)

#define Input_Output_i_8_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_8_512(out, v1);\
                                             Input_i_8_512(v1, in)

#define Input_Output_i_9_512(out,v1,in)      v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_9_512(out, v1);\
                                             Input_i_9_512(v1, in)

#define Input_Output_i_10_512(out,v1,in)     v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_10_512(out, v1);\
                                             Input_i_10_512(v1, in)

#define Input_Output_i_11_512(out,v1,in)     v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_11_512(out, v1);\
                                             Input_i_11_512(v1, in)

#define Input_Output_i_12_512(out,v1,in)     v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_12_512(out, v1);\
                                             Input_i_12_512(v1, in)

#define Input_Output_i_13_512(out,v1,in)     v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_13_512(out, v1);\
                                             Input_i_13_512(v1, in)

#define Input_Output_i_14_512(out,v1,in)     v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_14_512(out, v1);\
                                             Input_i_14_512(v1, in)

#define Input_Output_i_15_512(out,v1,in)     v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_15_512(out, v1);\
                                             Input_i_15_512(v1, in)

#define Input_Output_i_16_512(out,v1,in)     v1 = vrotate_high2lowi_512(v1);\
                                             Output_i_16_512(out, v1);\
                                             Input_i_16_512(v1, in)


static inline void transpose8x8_epi32_512_func(
    __m256i r0, __m256i r1, __m256i r2, __m256i r3,
    __m256i r4, __m256i r5, __m256i r6, __m256i r7,
    __m256i *o0, __m256i *o1, __m256i *o2, __m256i *o3,
    __m256i *o4, __m256i *o5, __m256i *o6, __m256i *o7)
{
    __m256i t0, t1, t2, t3, t4, t5, t6, t7;
    __m256i s0, s1, s2, s3, s4, s5, s6, s7;

    // Step 1: interleave adjacent rows by 32-bit elements.
    t0 = _mm256_unpacklo_epi32(r0, r1);
    t1 = _mm256_unpackhi_epi32(r0, r1);
    t2 = _mm256_unpacklo_epi32(r2, r3);
    t3 = _mm256_unpackhi_epi32(r2, r3);
    t4 = _mm256_unpacklo_epi32(r4, r5);
    t5 = _mm256_unpackhi_epi32(r4, r5);
    t6 = _mm256_unpacklo_epi32(r6, r7);
    t7 = _mm256_unpackhi_epi32(r6, r7);

    // Step 2: merge 32-bit pairs into 64-bit groups.
    s0 = _mm256_unpacklo_epi64(t0, t2);
    s1 = _mm256_unpackhi_epi64(t0, t2);
    s2 = _mm256_unpacklo_epi64(t1, t3);
    s3 = _mm256_unpackhi_epi64(t1, t3);
    s4 = _mm256_unpacklo_epi64(t4, t6);
    s5 = _mm256_unpackhi_epi64(t4, t6);
    s6 = _mm256_unpacklo_epi64(t5, t7);
    s7 = _mm256_unpackhi_epi64(t5, t7);

    // Step 3: exchange 128-bit halves to finish the 8x8 transpose.
    *o0 = _mm256_permute2x128_si256(s0, s4, 0x20);
    *o1 = _mm256_permute2x128_si256(s1, s5, 0x20);
    *o2 = _mm256_permute2x128_si256(s2, s6, 0x20);
    *o3 = _mm256_permute2x128_si256(s3, s7, 0x20);

    *o4 = _mm256_permute2x128_si256(s0, s4, 0x31);
    *o5 = _mm256_permute2x128_si256(s1, s5, 0x31);
    *o6 = _mm256_permute2x128_si256(s2, s6, 0x31);
    *o7 = _mm256_permute2x128_si256(s3, s7, 0x31);
}


static inline veci_512 combine_256_to_512_epi32(__m256i lo, __m256i hi)
{
    // Combine two 8-int vectors into one 16-int AVX-512 vector.
    return _mm512_inserti64x4(_mm512_castsi256_si512(lo), hi, 1);
}


static inline void transpose16x16_epi32_512_func(
    veci_512 r0,  veci_512 r1,  veci_512 r2,  veci_512 r3,
    veci_512 r4,  veci_512 r5,  veci_512 r6,  veci_512 r7,
    veci_512 r8,  veci_512 r9,  veci_512 r10, veci_512 r11,
    veci_512 r12, veci_512 r13, veci_512 r14, veci_512 r15,
    veci_512 *o0,  veci_512 *o1,  veci_512 *o2,  veci_512 *o3,
    veci_512 *o4,  veci_512 *o5,  veci_512 *o6,  veci_512 *o7,
    veci_512 *o8,  veci_512 *o9,  veci_512 *o10, veci_512 *o11,
    veci_512 *o12, veci_512 *o13, veci_512 *o14, veci_512 *o15)
{
    __m256i a0, a1, a2, a3, a4, a5, a6, a7;
    __m256i c0, c1, c2, c3, c4, c5, c6, c7;

    // Transpose rows 0..7, columns 0..7.
    transpose8x8_epi32_512_func(
        _mm512_castsi512_si256(r0),
        _mm512_castsi512_si256(r1),
        _mm512_castsi512_si256(r2),
        _mm512_castsi512_si256(r3),
        _mm512_castsi512_si256(r4),
        _mm512_castsi512_si256(r5),
        _mm512_castsi512_si256(r6),
        _mm512_castsi512_si256(r7),
        &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);

    // Transpose rows 8..15, columns 0..7.
    transpose8x8_epi32_512_func(
        _mm512_castsi512_si256(r8),
        _mm512_castsi512_si256(r9),
        _mm512_castsi512_si256(r10),
        _mm512_castsi512_si256(r11),
        _mm512_castsi512_si256(r12),
        _mm512_castsi512_si256(r13),
        _mm512_castsi512_si256(r14),
        _mm512_castsi512_si256(r15),
        &c0, &c1, &c2, &c3, &c4, &c5, &c6, &c7);

    // Build output columns 0..7 from the two lower-column 8x8 blocks.
    *o0 = combine_256_to_512_epi32(a0, c0);
    *o1 = combine_256_to_512_epi32(a1, c1);
    *o2 = combine_256_to_512_epi32(a2, c2);
    *o3 = combine_256_to_512_epi32(a3, c3);
    *o4 = combine_256_to_512_epi32(a4, c4);
    *o5 = combine_256_to_512_epi32(a5, c5);
    *o6 = combine_256_to_512_epi32(a6, c6);
    *o7 = combine_256_to_512_epi32(a7, c7);

    // Transpose rows 0..7, columns 8..15.
    transpose8x8_epi32_512_func(
        _mm512_extracti64x4_epi64(r0, 1),
        _mm512_extracti64x4_epi64(r1, 1),
        _mm512_extracti64x4_epi64(r2, 1),
        _mm512_extracti64x4_epi64(r3, 1),
        _mm512_extracti64x4_epi64(r4, 1),
        _mm512_extracti64x4_epi64(r5, 1),
        _mm512_extracti64x4_epi64(r6, 1),
        _mm512_extracti64x4_epi64(r7, 1),
        &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7);

    // Transpose rows 8..15, columns 8..15.
    transpose8x8_epi32_512_func(
        _mm512_extracti64x4_epi64(r8,  1),
        _mm512_extracti64x4_epi64(r9,  1),
        _mm512_extracti64x4_epi64(r10, 1),
        _mm512_extracti64x4_epi64(r11, 1),
        _mm512_extracti64x4_epi64(r12, 1),
        _mm512_extracti64x4_epi64(r13, 1),
        _mm512_extracti64x4_epi64(r14, 1),
        _mm512_extracti64x4_epi64(r15, 1),
        &c0, &c1, &c2, &c3, &c4, &c5, &c6, &c7);

    // Build output columns 8..15 from the two higher-column 8x8 blocks.
    *o8  = combine_256_to_512_epi32(a0, c0);
    *o9  = combine_256_to_512_epi32(a1, c1);
    *o10 = combine_256_to_512_epi32(a2, c2);
    *o11 = combine_256_to_512_epi32(a3, c3);
    *o12 = combine_256_to_512_epi32(a4, c4);
    *o13 = combine_256_to_512_epi32(a5, c5);
    *o14 = combine_256_to_512_epi32(a6, c6);
    *o15 = combine_256_to_512_epi32(a7, c7);
}


#define transposei_512(r0, r1, r2, r3, r4, r5, r6, r7, \
                       r8, r9, r10, r11, r12, r13, r14, r15, \
                       t0, t1, t2, t3, t4, t5, t6, t7, \
                       t8, t9, t10, t11, t12, t13, t14, t15) \
                       transpose16x16_epi32_512_func( \
                           r0, r1, r2, r3, r4, r5, r6, r7, \
                           r8, r9, r10, r11, r12, r13, r14, r15, \
                           &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, \
                           &t8, &t9, &t10, &t11, &t12, &t13, &t14, &t15)

void vectime_512(int* A, int NX, int NY, int T);
void vectime_extra_array_512(int* A, int NX, int NY, int T);
