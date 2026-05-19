#include "define_512.h"

void vectime_transpose_boundary_extra_array_512_hybrid(double* A, int NX, int NY, int NZ, int T) {
	if (back_scalar_512(NX, NY, NZ, T)) {
		naive_scalar(A, NX, NY, NZ, T);
		return;
	}
    
    double (* B)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;

    long int i, j, k;

    double tmp[VVECLEN];
    int tt, t = 0, x, xx, y, yy, z, zz;

    vec_512 v_x_plus_1, v_x_minus_1;
    vec_512 v_y_plus_1, v_y_minus_1;
    vec_512 v_center_0, v_center_1, v_center_2, v_center_3, v_center_4;
    vec_512 v_center_5, v_center_6, v_center_7, v_center_8, v_center_9; 

    vec_512 in, out;
    SET_COFF_512;
    
    int Y_EXT = NY + 2 * YSTART;
    double (* AV) [Y_EXT][NZ][VVECLEN] = (double(*)[Y_EXT][NZ][VVECLEN])alloc_extra_array_512(sizeof(double) * Y_EXT * NZ * VVECLEN * 3);

    double (* BV0) [NZ][VVECLEN] = (double (*) [NZ][VVECLEN]) AV;
    double (* BV1) [NZ][VVECLEN] = (double (*) [NZ][VVECLEN]) (AV + 1);
    double (* BV2) [NZ][VVECLEN] = (double (*) [NZ][VVECLEN]) (AV + 2);

    double (* Btmp [3]) [NZ][VVECLEN]  = {BV0, BV1, BV2};

    for ( tt = 0; tt <= T - VVECLEN; tt += VVECLEN){  
        for( t = tt ; t < tt + VVECLEN - 1 ; t++){      //head
            for ( x = XSTART; x < XSTART + STRIDE * (VVECLEN - 1 - (t - tt)); x++) {
                for ( y = YSTART; y < NY + YSTART; y++) {
                    #pragma ivdep
                    #pragma vector always
                    // AVX-256 7-point Z-SIMD: z-1+z+1+y-1+y+1+x-1+x+1 order
				__m256d vc0 = _mm256_set1_pd(C0);
				__m256d vc1 = _mm256_set1_pd(C1);
				for (z = ZSTART; z + 3 < NZ + ZSTART; z += 4){
					__m256d vm  = _mm256_loadu_pd(&B[t%2][x][y][z]);
					__m256d vzm = _mm256_loadu_pd(&B[t%2][x][y][z-1]);
					__m256d vzp = _mm256_loadu_pd(&B[t%2][x][y][z+1]);
					__m256d vym = _mm256_loadu_pd(&B[t%2][x][y-1][z]);
					__m256d vyp = _mm256_loadu_pd(&B[t%2][x][y+1][z]);
					__m256d vxm = _mm256_loadu_pd(&B[t%2][x-1][y][z]);
					__m256d vxp = _mm256_loadu_pd(&B[t%2][x+1][y][z]);
					__m256d stub = _mm256_add_pd(vzm, vzp);
					stub = _mm256_add_pd(stub, vym);
					stub = _mm256_add_pd(stub, vyp);
					stub = _mm256_add_pd(stub, vxm);
					stub = _mm256_add_pd(stub, vxp);
					__m256d res = _mm256_fmadd_pd(vc0, vm, _mm256_mul_pd(vc1, stub));
					_mm256_storeu_pd(&B[(t+1)%2][x][y][z], res);
				}
				for (; z < NZ + ZSTART; z++){
					Compute_scalar(B, t, x, y, z);
				}   
                }       
            }
        }
        t = tt;

        //初始转置
        for(x = XSTART - XSLOPE; x <= XSTART + XSLOPE; x++){
            for ( y = YSTART - YSLOPE; y < NY + YSTART + YSLOPE ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN ) {
                    vload_512(v_center_0, B[t%2    ][x + STRIDE *7][y][z]);
                    vload_512(v_center_1, B[(t+1)%2][x + STRIDE *6][y][z]);
                    vload_512(v_center_2, B[t%2    ][x + STRIDE *5][y][z]);
                    vload_512(v_center_3, B[(t+1)%2][x + STRIDE *4][y][z]);
                    vload_512(v_center_4, B[t%2    ][x + STRIDE *3][y][z]);
                    vload_512(v_center_5, B[(t+1)%2][x + STRIDE *2][y][z]);
                    vload_512(v_center_6, B[t%2    ][x + STRIDE *1][y][z]);
                    vload_512(v_center_7, B[(t+1)%2][x            ][y][z]);
                    
                    transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);
                    
                    int x_idx = x - XSTART + XSLOPE;
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 0][0], v_center_0);
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 1][0], v_center_1);
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 2][0], v_center_2);
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 3][0], v_center_3);
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 4][0], v_center_4);
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 5][0], v_center_5);
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 6][0], v_center_6);
                    vstore_512(Btmp[x_idx][y][z - ZSTART + 7][0], v_center_7);
                }
            }
        } 
        
        // 双面转置
        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VVECLEN ; x ++){

            int x_new = x + STRIDE; 

            // [下边界 y = YSTART - 1]
            y = YSTART - YSLOPE;
            for ( z = ZSTART ; z <= NZ + ZSTART - VVECLEN; z += VVECLEN){ 
                vload_512(v_center_0, B[t%2    ][x_new + STRIDE *7][y][z]);
                vload_512(v_center_1, B[(t+1)%2][x_new + STRIDE *6][y][z]);
                vload_512(v_center_2, B[t%2    ][x_new + STRIDE *5][y][z]);
                vload_512(v_center_3, B[(t+1)%2][x_new + STRIDE *4][y][z]);
                vload_512(v_center_4, B[t%2    ][x_new + STRIDE *3][y][z]);
                vload_512(v_center_5, B[(t+1)%2][x_new + STRIDE *2][y][z]);
                vload_512(v_center_6, B[t%2    ][x_new + STRIDE *1][y][z]);
                vload_512(v_center_7, B[(t+1)%2][x_new            ][y][z]);
                
                transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);
                
                vstore_512(BV0[y][z - ZSTART + 0][0], v_center_0);
                vstore_512(BV0[y][z - ZSTART + 1][0], v_center_1);
                vstore_512(BV0[y][z - ZSTART + 2][0], v_center_2);
                vstore_512(BV0[y][z - ZSTART + 3][0], v_center_3);
                vstore_512(BV0[y][z - ZSTART + 4][0], v_center_4);
                vstore_512(BV0[y][z - ZSTART + 5][0], v_center_5);
                vstore_512(BV0[y][z - ZSTART + 6][0], v_center_6);
                vstore_512(BV0[y][z - ZSTART + 7][0], v_center_7);
            }

            // [上边界 y = NY + YSTART]
            y = NY + YSTART;
            for ( z = ZSTART ; z <= NZ + ZSTART - VVECLEN; z += VVECLEN){ 
                vload_512(v_center_0, B[t%2    ][x_new + STRIDE *7][y][z]);
                vload_512(v_center_1, B[(t+1)%2][x_new + STRIDE *6][y][z]);
                vload_512(v_center_2, B[t%2    ][x_new + STRIDE *5][y][z]);
                vload_512(v_center_3, B[(t+1)%2][x_new + STRIDE *4][y][z]);
                vload_512(v_center_4, B[t%2    ][x_new + STRIDE *3][y][z]);
                vload_512(v_center_5, B[(t+1)%2][x_new + STRIDE *2][y][z]);
                vload_512(v_center_6, B[t%2    ][x_new + STRIDE *1][y][z]);
                vload_512(v_center_7, B[(t+1)%2][x_new            ][y][z]);
                
                transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);
                
                vstore_512(BV0[y][z - ZSTART + 0][0], v_center_0);
                vstore_512(BV0[y][z - ZSTART + 1][0], v_center_1);
                vstore_512(BV0[y][z - ZSTART + 2][0], v_center_2);
                vstore_512(BV0[y][z - ZSTART + 3][0], v_center_3);
                vstore_512(BV0[y][z - ZSTART + 4][0], v_center_4);
                vstore_512(BV0[y][z - ZSTART + 5][0], v_center_5);
                vstore_512(BV0[y][z - ZSTART + 6][0], v_center_6);
                vstore_512(BV0[y][z - ZSTART + 7][0], v_center_7);
            }

            // 核心计算区
            for ( y = YSTART; y < NY + YSTART; y++ ) {
                z = ZSTART;
                
                vset_3d_512(v_center_0, x, y, z-ZSLOPE); 
                vload_512(v_center_1, BV1[y][0][0]); 

                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN - VVECLEN; z += VVECLEN) {
                
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]); 

                    // 1st
                    vload_512(v_center_2 , BV1[y][z - ZSTART + 1][0]);                                                                        
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART][0]));
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART][0])));   
                    Input_Output_1_512(out, v_center_0, in); vstore_512(BV0[y][z - ZSTART][0], v_center_0);
                    
                    // 2nd
                    vload_512(v_center_3 , BV1[y][z - ZSTART + 2][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 1][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 1][0]));
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 1][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 1][0])));
                    Input_Output_2_512(out, v_center_1, in); vstore_512(BV0[y][z - ZSTART + 1][0], v_center_1); 
                    
                    // 3rd
                    vload_512(v_center_4 , BV1[y][z - ZSTART + 3][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 2][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 2][0]));
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 2][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 2][0]))); 
                    Input_Output_3_512(out, v_center_2, in); vstore_512(BV0[y][z - ZSTART + 2][0], v_center_2); 
                    
                    // 4th
                    vload_512(v_center_5 , BV1[y][z - ZSTART + 4][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 3][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 3][0]));
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 3][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 3][0]))); 
                    Input_Output_4_512(out, v_center_3, in); vstore_512(BV0[y][z - ZSTART + 3][0], v_center_3); 
                    
                    // 5th
                    vload_512(v_center_6 , BV1[y][z - ZSTART + 5][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 4][0]));    
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 4][0]));    
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 4][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 4][0])));
                    Input_Output_5_512(out, v_center_4, in); vstore_512(BV0[y][z - ZSTART + 4][0], v_center_4); 

                    // 6th
                    vload_512(v_center_7 , BV1[y][z - ZSTART + 6][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 5][0]));    
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 5][0]));    
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 5][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 5][0]))); 
                    Input_Output_6_512(out, v_center_5, in); vstore_512(BV0[y][z - ZSTART + 5][0], v_center_5); 

                    // 7th
                    vload_512(v_center_8 , BV1[y][z - ZSTART + 7][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 6][0]));    
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 6][0]));    
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 6][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 6][0]))); 
                    Input_Output_7_512(out, v_center_6, in); vstore_512(BV0[y][z - ZSTART + 6][0], v_center_6); 

                    // 8th
                    v_center_9 = _mm512_loadu_pd(&(BV1[y][z - ZSTART + 8][0])); 
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 7][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 7][0]));
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 7][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 7][0]))); 
                    Input_Output_8_512(out, v_center_7, in); vstore_512(BV0[y][z - ZSTART + 7][0], v_center_7); 
                    
                    vstore_512(B[(t)%2][x][y][z], out);
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                }

                // 最后一个 Z 向量越界保护，外层读取主存，补齐 z += VVECLEN
                if ( z <= NZ + ZSTART - VVECLEN ) {
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]);   

                    // 1st
                    vload_512(v_center_2 , BV1[y][z - ZSTART + 1][0]);                                                                        
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART][0]));
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART][0])));   
                    Input_Output_1_512(out, v_center_0, in); vstore_512(BV0[y][z - ZSTART][0], v_center_0);
                    
                    // 2nd
                    vload_512(v_center_3 , BV1[y][z - ZSTART + 2][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 1][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 1][0]));
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 1][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 1][0])));
                    Input_Output_2_512(out, v_center_1, in); vstore_512(BV0[y][z - ZSTART + 1][0], v_center_1); 
                    
                    // 3rd
                    vload_512(v_center_4 , BV1[y][z - ZSTART + 3][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 2][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 2][0]));
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 2][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 2][0]))); 
                    Input_Output_3_512(out, v_center_2, in); vstore_512(BV0[y][z - ZSTART + 2][0], v_center_2); 
                    
                    // 4th
                    vload_512(v_center_5 , BV1[y][z - ZSTART + 4][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 3][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 3][0]));
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 3][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 3][0]))); 
                    Input_Output_4_512(out, v_center_3, in); vstore_512(BV0[y][z - ZSTART + 3][0], v_center_3); 
                    
                    // 5th
                    vload_512(v_center_6 , BV1[y][z - ZSTART + 5][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 4][0]));    
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 4][0]));    
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 4][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 4][0])));
                    Input_Output_5_512(out, v_center_4, in); vstore_512(BV0[y][z - ZSTART + 4][0], v_center_4); 

                    // 6th
                    vload_512(v_center_7 , BV1[y][z - ZSTART + 6][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 5][0]));    
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 5][0]));    
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 5][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 5][0]))); 
                    Input_Output_6_512(out, v_center_5, in); vstore_512(BV0[y][z - ZSTART + 5][0], v_center_5); 

                    // 7th
                    vload_512(v_center_8 , BV1[y][z - ZSTART + 7][0]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 6][0]));    
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 6][0]));    
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 6][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 6][0]))); 
                    Input_Output_7_512(out, v_center_6, in); vstore_512(BV0[y][z - ZSTART + 6][0], v_center_6); 

                    // 8th
                    v_center_9 = ( z > NZ + ZSTART - VVECLEN - VVECLEN ) ? setv_3d_512(x, y, z + 8) : _mm512_loadu_pd(&(BV1[y][z - ZSTART + 8][0]));
                    v_y_minus_1 = _mm512_loadu_pd(&(BV1[y - YSLOPE][z - ZSTART + 7][0]));
                    v_y_plus_1  = _mm512_loadu_pd(&(BV1[y + YSLOPE][z - ZSTART + 7][0]));
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, v_y_plus_1, _mm512_loadu_pd(&(BV0[y][z - ZSTART + 7][0])), _mm512_loadu_pd(&(BV2[y][z - ZSTART + 7][0]))); 
                    Input_Output_8_512(out, v_center_7, in); vstore_512(BV0[y][z - ZSTART + 7][0], v_center_7); 
                    
                    vstore_512(B[(t)%2][x][y][z], out);
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;

                    z += VVECLEN;
                }

                // 标量收尾 
                for ( ; z < NZ + ZSTART; z++) {
                    vset_3d_512(v_center_2 , x, y, z + 1);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, setv_3d_512(x, y - 1, z), setv_3d_512(x, y + 1, z), setv_3d_512(x - 1, y, z), setv_3d_512(x + 1, y, z));
                    
                    _mm512_storeu_pd(tmp, v_center_0);
                    B[(t+1)%2][x+STRIDE*7][y][z] = tmp[0]; 
                    B[t%2    ][x+STRIDE*6][y][z] = tmp[1]; 
                    B[(t+1)%2][x+STRIDE*5][y][z] = tmp[2]; 
                    B[t%2    ][x+STRIDE*4][y][z] = tmp[3]; 
                    B[(t+1)%2][x+STRIDE*3][y][z] = tmp[4]; 
                    B[t%2    ][x+STRIDE*2][y][z] = tmp[5]; 
                    B[(t+1)%2][x+STRIDE*1][y][z] = tmp[6]; 
                    B[t%2    ][x         ][y][z] = tmp[7];

                    v_center_0 = v_center_1; 
                    v_center_1 = v_center_2;
                }
            }

            // 滚动切片轮替
            Btmp[0] = BV0;
            BV0 = BV1;
            BV1 = BV2;
            BV2 = Btmp[0];     
        }
        


        Btmp [0] = BV0;
        Btmp [1] = BV1;
        Btmp [2] = BV2;
        
        for(; x < NX + XSTART - STRIDE * VVECLEN + 1 + 3; x++){           
            for ( y = YSTART ; y < NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                    vload_512(v_center_0, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 0][0]);
                    vload_512(v_center_1, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 1][0]);
                    vload_512(v_center_2, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 2][0]);
                    vload_512(v_center_3, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 3][0]);    
                    vload_512(v_center_4, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 4][0]);
                    vload_512(v_center_5, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 5][0]);
                    vload_512(v_center_6, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 6][0]);
                    vload_512(v_center_7, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][y][z - ZSTART + 7][0]);
                    
                    transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);
                    
                    vstore_512( B[t%2    ][x - XSLOPE + STRIDE *7][y][z], v_center_0);
                    vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *6][y][z], v_center_1);
                    vstore_512( B[t%2    ][x - XSLOPE + STRIDE *5][y][z], v_center_2);
                    vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *4][y][z], v_center_3);
                    vstore_512( B[t%2    ][x - XSLOPE + STRIDE *3][y][z], v_center_4);
                    vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *2][y][z], v_center_5);
                    vstore_512( B[t%2    ][x - XSLOPE + STRIDE *1][y][z], v_center_6);
                    vstore_512( B[(t+1)%2][x - XSLOPE            ][y][z], v_center_7);
                }
            }
        }            

        // 纯标量扫尾区
        xx = NX + XSTART - STRIDE * VVECLEN + 1;
        for( t = tt ; t < tt + VVECLEN ; t++){   
            for ( x = xx + STRIDE * (VVECLEN - 1 - (t - tt)); x < NX + XSTART; x++) {
                for ( y = YSTART; y < NY + YSTART; y++) {
                    #pragma ivdep
                    #pragma vector always
                    // AVX-256 7-point Z-SIMD: z-1+z+1+y-1+y+1+x-1+x+1 order
				__m256d vc0 = _mm256_set1_pd(C0);
				__m256d vc1 = _mm256_set1_pd(C1);
				for (z = ZSTART; z + 3 < NZ + ZSTART; z += 4){
					__m256d vm  = _mm256_loadu_pd(&B[t%2][x][y][z]);
					__m256d vzm = _mm256_loadu_pd(&B[t%2][x][y][z-1]);
					__m256d vzp = _mm256_loadu_pd(&B[t%2][x][y][z+1]);
					__m256d vym = _mm256_loadu_pd(&B[t%2][x][y-1][z]);
					__m256d vyp = _mm256_loadu_pd(&B[t%2][x][y+1][z]);
					__m256d vxm = _mm256_loadu_pd(&B[t%2][x-1][y][z]);
					__m256d vxp = _mm256_loadu_pd(&B[t%2][x+1][y][z]);
					__m256d stub = _mm256_add_pd(vzm, vzp);
					stub = _mm256_add_pd(stub, vym);
					stub = _mm256_add_pd(stub, vyp);
					stub = _mm256_add_pd(stub, vxm);
					stub = _mm256_add_pd(stub, vxp);
					__m256d res = _mm256_fmadd_pd(vc0, vm, _mm256_mul_pd(vc1, stub));
					_mm256_storeu_pd(&B[(t+1)%2][x][y][z], res);
				}
				for (; z < NZ + ZSTART; z++){
					Compute_scalar(B, t, x, y, z);
				}   
                }   
            }
        }
    }
    
    // Extra points
    for ( ; t < T; t++){
        for (x = XSTART; x < NX + XSTART; x++) {
            for ( y = YSTART; y < NY + YSTART; y++) {
                #pragma ivdep
                #pragma vector always
                // AVX-256 7-point Z-SIMD: z-1+z+1+y-1+y+1+x-1+x+1 order
				__m256d vc0 = _mm256_set1_pd(C0);
				__m256d vc1 = _mm256_set1_pd(C1);
				for (z = ZSTART; z + 3 < NZ + ZSTART; z += 4){
					__m256d vm  = _mm256_loadu_pd(&B[t%2][x][y][z]);
					__m256d vzm = _mm256_loadu_pd(&B[t%2][x][y][z-1]);
					__m256d vzp = _mm256_loadu_pd(&B[t%2][x][y][z+1]);
					__m256d vym = _mm256_loadu_pd(&B[t%2][x][y-1][z]);
					__m256d vyp = _mm256_loadu_pd(&B[t%2][x][y+1][z]);
					__m256d vxm = _mm256_loadu_pd(&B[t%2][x-1][y][z]);
					__m256d vxp = _mm256_loadu_pd(&B[t%2][x+1][y][z]);
					__m256d stub = _mm256_add_pd(vzm, vzp);
					stub = _mm256_add_pd(stub, vym);
					stub = _mm256_add_pd(stub, vyp);
					stub = _mm256_add_pd(stub, vxm);
					stub = _mm256_add_pd(stub, vxp);
					__m256d res = _mm256_fmadd_pd(vc0, vm, _mm256_mul_pd(vc1, stub));
					_mm256_storeu_pd(&B[(t+1)%2][x][y][z], res);
				}
				for (; z < NZ + ZSTART; z++){
					Compute_scalar(B, t, x, y, z);
				}
            }
        }
    }   
    
    free_extra_array(AV);
}
