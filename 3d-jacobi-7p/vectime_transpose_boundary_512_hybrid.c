#include "define_512.h"

void vectime_transpose_boundary_512_hybrid(double* A, int NX, int NY, int NZ, int T) {
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


        for(x = XSTART - XSLOPE; x <= XSTART + XSLOPE; x++){
            for ( y = YSTART - YSLOPE; y <= NY + YSTART ; y ++ ){
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
                    
                    vstore_512( B[(t+1)%2][x            ][y][z], v_center_0);
                    vstore_512( B[t%2    ][x + STRIDE *1][y][z], v_center_1);
                    vstore_512( B[(t+1)%2][x + STRIDE *2][y][z], v_center_2);
                    vstore_512( B[t%2    ][x + STRIDE *3][y][z], v_center_3);
                    vstore_512( B[(t+1)%2][x + STRIDE *4][y][z], v_center_4);
                    vstore_512( B[t%2    ][x + STRIDE *5][y][z], v_center_5);
                    vstore_512( B[(t+1)%2][x + STRIDE *6][y][z], v_center_6);
                    vstore_512( B[t%2    ][x + STRIDE *7][y][z], v_center_7);
                }
            }
        }

        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VVECLEN ; x ++){
            y = YSTART - YSLOPE;
            for ( z = ZSTART ; z <= NZ + ZSTART - VVECLEN; z += VVECLEN){ 
                
                vload_512(v_center_0, B[t%2    ][x + XSLOPE + 1 + STRIDE *7][y][z]);
                vload_512(v_center_1, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *6][y][z]);
                vload_512(v_center_2, B[t%2    ][x + XSLOPE + 1 + STRIDE *5][y][z]);
                vload_512(v_center_3, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *4][y][z]);
                vload_512(v_center_4, B[t%2    ][x + XSLOPE + 1 + STRIDE *3][y][z]);
                vload_512(v_center_5, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][y][z]);
                vload_512(v_center_6, B[t%2    ][x + XSLOPE + 1 + STRIDE *1][y][z]);
                vload_512(v_center_7, B[(t+1)%2][x + XSLOPE + 1            ][y][z]);

                transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);

                vstore_512( B[(t+1)%2][x + XSLOPE + 1            ][y][z], v_center_0);
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *1][y][z], v_center_1);
                vstore_512( B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][y][z], v_center_2);
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *3][y][z], v_center_3);
                vstore_512( B[(t+1)%2][x + XSLOPE + 1 + STRIDE *4][y][z], v_center_4);
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *5][y][z], v_center_5);
                vstore_512( B[(t+1)%2][x + XSLOPE + 1 + STRIDE *6][y][z], v_center_6);
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *7][y][z], v_center_7);

                vload_512(v_center_0, B[t%2    ][x + XSLOPE + 1 + STRIDE *7][y + NY + YSLOPE][z]); 
                vload_512(v_center_1, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *6][y + NY + YSLOPE][z]);
                vload_512(v_center_2, B[t%2    ][x + XSLOPE + 1 + STRIDE *5][y + NY + YSLOPE][z]);
                vload_512(v_center_3, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *4][y + NY + YSLOPE][z]);
                vload_512(v_center_4, B[t%2    ][x + XSLOPE + 1 + STRIDE *3][y + NY + YSLOPE][z]);
                vload_512(v_center_5, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][y + NY + YSLOPE][z]);
                vload_512(v_center_6, B[t%2    ][x + XSLOPE + 1 + STRIDE *1][y + NY + YSLOPE][z]);
                vload_512(v_center_7, B[(t+1)%2][x + XSLOPE + 1            ][y + NY + YSLOPE][z]);

                transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);

                vstore_512( B[(t+1)%2][x + XSLOPE + 1            ][y + NY + YSLOPE][z], v_center_0 );
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *1][y + NY + YSLOPE][z], v_center_1 );
                vstore_512( B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][y + NY + YSLOPE][z], v_center_2 );
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *3][y + NY + YSLOPE][z], v_center_3 );
                vstore_512( B[(t+1)%2][x + XSLOPE + 1 + STRIDE *4][y + NY + YSLOPE][z], v_center_4 );
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *5][y + NY + YSLOPE][z], v_center_5 );
                vstore_512( B[(t+1)%2][x + XSLOPE + 1 + STRIDE *6][y + NY + YSLOPE][z], v_center_6 );
                vstore_512( B[t%2    ][x + XSLOPE + 1 + STRIDE *7][y + NY + YSLOPE][z], v_center_7 );
            }

            // CORE Y LOOP
            for ( y = YSTART; y < NY + YSTART; y++ ) {
                z = ZSTART;
                vloadset_512(v_center_0, 
                    B[(t+1)%2][x         ][y][z-ZSLOPE], 
                    B[t%2    ][x+STRIDE*1][y][z-ZSLOPE], 
                    B[(t+1)%2][x+STRIDE*2][y][z-ZSLOPE], 
                    B[t%2    ][x+STRIDE*3][y][z-ZSLOPE], 
                    B[(t+1)%2][x+STRIDE*4][y][z-ZSLOPE], 
                    B[t%2    ][x+STRIDE*5][y][z-ZSLOPE], 
                    B[(t+1)%2][x+STRIDE*6][y][z-ZSLOPE], 
                    B[t%2    ][x+STRIDE*7][y][z-ZSLOPE]);
                                    
                vload_512(v_center_1, B[(t+1)%2  ][x][y][z]);   

                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]); 

                    vload_512(v_center_2 , B[(t)%2  ][x+STRIDE][y][z]);                                                                      
                    v_y_minus_1 = _mm512_loadu_pd(&B[(t+1)%2][x][y-YSLOPE][z]);
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t+1)%2][x][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE][y][z]);
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE][y][z]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); 
                    Input_Output_1_512(out, v_center_0, in);
                    vstore_512(B[(t+1)%2][x+STRIDE][y][z], v_center_0); 
                    
                    vload_512(v_center_3 , B[(t+1)%2][x+STRIDE*2][y][z]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE][y-YSLOPE][z]);
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE][y][z]); 
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE][y][z]);
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); 
                    Input_Output_2_512(out, v_center_1, in);
                    vstore_512(B[(t)%2][x+STRIDE*2][y][z], v_center_1); 
                    
                    vload_512(v_center_4 , B[(t)%2  ][x+STRIDE*3][y][z]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z]);
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*2][y][z]); 
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*2][y][z]);
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); 
                    Input_Output_3_512(out, v_center_2, in);    
                    vstore_512(B[(t+1)%2  ][x+STRIDE*3][y][z], v_center_2); 
                    
                    vload_512(v_center_5 , B[(t+1)%2][x+STRIDE*4][y][z]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE*3][y-YSLOPE][z]);
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE*3][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*3][y][z]); 
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*3][y][z]);
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); 
                    Input_Output_4_512(out, v_center_3, in);    
                    vstore_512(B[(t)%2][x+STRIDE*4][y][z], v_center_3); 
                    
                    vload_512(v_center_6 , B[(t)%2  ][x+STRIDE*5][y][z]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z]);    
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*4][y][z]); 
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*4][y][z]);
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_5_512(out, v_center_4, in);    
                    vstore_512(B[(t+1)%2][x+STRIDE*5][y][z], v_center_4); 

                    vload_512(v_center_7 , B[(t+1)%2][x+STRIDE*6][y][z]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE*5][y-YSLOPE][z]);    
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE*5][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*5][y][z]); 
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*5][y][z]);
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); 
                    Input_Output_6_512(out, v_center_5, in);    
                    vstore_512(B[(t)%2][x+STRIDE*6][y][z], v_center_5); 

                    vload_512(v_center_8 , B[(t)%2  ][x+STRIDE*7][y][z]);                  
                    v_y_minus_1 = _mm512_loadu_pd(&B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z]);    
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*6][y][z]); 
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*6][y][z]);
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); 
                    Input_Output_7_512(out, v_center_6, in);    
                    vstore_512(B[(t+1)%2][x+STRIDE*7][y][z], v_center_6); 

                    v_center_9 = ( z > NZ + ZSTART - VVECLEN - VVECLEN ) ? 
                        _mm512_set_pd(
                            B[(t+1)%2][x         ][y][z+VVECLEN], 
                            B[t%2    ][x+STRIDE*1][y][z+VVECLEN], 
                            B[(t+1)%2][x+STRIDE*2][y][z+VVECLEN], 
                            B[t%2    ][x+STRIDE*3][y][z+VVECLEN], 
                            B[(t+1)%2][x+STRIDE*4][y][z+VVECLEN], 
                            B[t%2    ][x+STRIDE*5][y][z+VVECLEN], 
                            B[(t+1)%2][x+STRIDE*6][y][z+VVECLEN], 
                            B[t%2    ][x+STRIDE*7][y][z+VVECLEN]
                        ) : _mm512_loadu_pd(&B[(t+1)%2][x][y][z+VVECLEN]);

                    v_y_minus_1 = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE*7][y-YSLOPE][z]); 
                    v_y_plus_1  = _mm512_loadu_pd(&B[(t  )%2][x+STRIDE*7][y+YSLOPE][z]); 
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*7][y][z]); 
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*7][y][z]);
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); 
                    Input_Output_8_512(out, v_center_7, in);    
                    vstore_512(B[(t)%2][x+STRIDE*8][y][z], v_center_7); 
                    
                    vstore_512(B[(t)%2][x][y][z], out);
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                }

                // 标量扫尾
                for ( ; z < NZ + ZSTART; z++) {
                    vloadset_512(v_y_minus_1, 
                        B[(t+1)%2][x         ][y-1][z], B[t%2    ][x+STRIDE*1][y-1][z], 
                        B[(t+1)%2][x+STRIDE*2][y-1][z], B[t%2    ][x+STRIDE*3][y-1][z], 
                        B[(t+1)%2][x+STRIDE*4][y-1][z], B[t%2    ][x+STRIDE*5][y-1][z], 
                        B[(t+1)%2][x+STRIDE*6][y-1][z], B[t%2    ][x+STRIDE*7][y-1][z]);
                    vloadset_512(v_y_plus_1, 
                        B[(t+1)%2][x         ][y+1][z], B[t%2    ][x+STRIDE*1][y+1][z], 
                        B[(t+1)%2][x+STRIDE*2][y+1][z], B[t%2    ][x+STRIDE*3][y+1][z], 
                        B[(t+1)%2][x+STRIDE*4][y+1][z], B[t%2    ][x+STRIDE*5][y+1][z], 
                        B[(t+1)%2][x+STRIDE*6][y+1][z], B[t%2    ][x+STRIDE*7][y+1][z]);
                    vloadset_512(v_x_minus_1, 
                        B[(t+1)%2][x-1         ][y][z], B[t%2    ][x-1+STRIDE*1][y][z], 
                        B[(t+1)%2][x-1+STRIDE*2][y][z], B[t%2    ][x-1+STRIDE*3][y][z], 
                        B[(t+1)%2][x-1+STRIDE*4][y][z], B[t%2    ][x-1+STRIDE*5][y][z], 
                        B[(t+1)%2][x-1+STRIDE*6][y][z], B[t%2    ][x-1+STRIDE*7][y][z]);
                    vloadset_512(v_x_plus_1, 
                        B[(t+1)%2][x+1         ][y][z], B[t%2    ][x+1+STRIDE*1][y][z], 
                        B[(t+1)%2][x+1+STRIDE*2][y][z], B[t%2    ][x+1+STRIDE*3][y][z], 
                        B[(t+1)%2][x+1+STRIDE*4][y][z], B[t%2    ][x+1+STRIDE*5][y][z], 
                        B[(t+1)%2][x+1+STRIDE*6][y][z], B[t%2    ][x+1+STRIDE*7][y][z]);
                    vloadset_512(v_center_2, 
                        B[(t+1)%2][x         ][y][z+1], B[t%2    ][x+STRIDE*1][y][z+1], 
                        B[(t+1)%2][x+STRIDE*2][y][z+1], B[t%2    ][x+STRIDE*3][y][z+1], 
                        B[(t+1)%2][x+STRIDE*4][y][z+1], B[t%2    ][x+STRIDE*5][y][z+1], 
                        B[(t+1)%2][x+STRIDE*6][y][z+1], B[t%2    ][x+STRIDE*7][y][z+1]);

                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    
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

            y = YSTART - YSLOPE;
            for ( z = ZSTART ; z <= NZ + ZSTART - VVECLEN; z += VVECLEN){ 

                vload_512(v_center_0, B[(t+1)%2][x - XSLOPE            ][y][z]);
                vload_512(v_center_1, B[t%2    ][x - XSLOPE + STRIDE *1][y][z]);
                vload_512(v_center_2, B[(t+1)%2][x - XSLOPE + STRIDE *2][y][z]);
                vload_512(v_center_3, B[t%2    ][x - XSLOPE + STRIDE *3][y][z]);
                vload_512(v_center_4, B[(t+1)%2][x - XSLOPE + STRIDE *4][y][z]);
                vload_512(v_center_5, B[t%2    ][x - XSLOPE + STRIDE *5][y][z]);
                vload_512(v_center_6, B[(t+1)%2][x - XSLOPE + STRIDE *6][y][z]);
                vload_512(v_center_7, B[t%2    ][x - XSLOPE + STRIDE *7][y][z]);

                transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);

                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *7][y][z], v_center_0);
                vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *6][y][z], v_center_1);
                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *5][y][z], v_center_2);
                vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *4][y][z], v_center_3);
                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *3][y][z], v_center_4);
                vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *2][y][z], v_center_5);
                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *1][y][z], v_center_6);
                vstore_512( B[(t+1)%2][x - XSLOPE            ][y][z], v_center_7);

                vload_512(v_center_0, B[(t+1)%2][x - XSLOPE            ][y + NY + YSLOPE][z]);
                vload_512(v_center_1, B[t%2    ][x - XSLOPE + STRIDE *1][y + NY + YSLOPE][z]);
                vload_512(v_center_2, B[(t+1)%2][x - XSLOPE + STRIDE *2][y + NY + YSLOPE][z]);
                vload_512(v_center_3, B[t%2    ][x - XSLOPE + STRIDE *3][y + NY + YSLOPE][z]);
                vload_512(v_center_4, B[(t+1)%2][x - XSLOPE + STRIDE *4][y + NY + YSLOPE][z]);
                vload_512(v_center_5, B[t%2    ][x - XSLOPE + STRIDE *5][y + NY + YSLOPE][z]);
                vload_512(v_center_6, B[(t+1)%2][x - XSLOPE + STRIDE *6][y + NY + YSLOPE][z]);
                vload_512(v_center_7, B[t%2    ][x - XSLOPE + STRIDE *7][y + NY + YSLOPE][z]);

                transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);

                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *7][y + NY + YSLOPE][z], v_center_0 );
                vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *6][y + NY + YSLOPE][z], v_center_1 );
                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *5][y + NY + YSLOPE][z], v_center_2 );
                vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *4][y + NY + YSLOPE][z], v_center_3 );
                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *3][y + NY + YSLOPE][z], v_center_4 );
                vstore_512( B[(t+1)%2][x - XSLOPE + STRIDE *2][y + NY + YSLOPE][z], v_center_5 );
                vstore_512( B[t%2    ][x - XSLOPE + STRIDE *1][y + NY + YSLOPE][z], v_center_6 );
                vstore_512( B[(t+1)%2][x - XSLOPE            ][y + NY + YSLOPE][z], v_center_7 );
            }
        }
        

        // Tail Transpose
        for(; x < NX + XSTART - STRIDE * VVECLEN + 1 + 3; x++){          
            for ( y = YSTART - YSLOPE ; y <= NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                    vload_512(v_center_0, B[(t+1)%2][x - XSLOPE            ][y][z]);
                    vload_512(v_center_1, B[t%2    ][x - XSLOPE + STRIDE *1][y][z]);
                    vload_512(v_center_2, B[(t+1)%2][x - XSLOPE + STRIDE *2][y][z]);
                    vload_512(v_center_3, B[t%2    ][x - XSLOPE + STRIDE *3][y][z]);
                    vload_512(v_center_4, B[(t+1)%2][x - XSLOPE + STRIDE *4][y][z]);
                    vload_512(v_center_5, B[t%2    ][x - XSLOPE + STRIDE *5][y][z]);
                    vload_512(v_center_6, B[(t+1)%2][x - XSLOPE + STRIDE *6][y][z]);
                    vload_512(v_center_7, B[t%2    ][x - XSLOPE + STRIDE *7][y][z]);
                    
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
        
        // Tail Computations
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
        for ( y = 0; y < NY + YSTART * 2; y++) {
            for ( z = 0; z < NZ + ZSTART * 2 ; z++) {
                B[1][XSTART - XSLOPE][y][z] = B[0][XSTART - XSLOPE][y][z]; 
            }
        }
    }
    
    // Extra points
    for ( ; t < T; t++){
        for (x = XSTART; x < NX + XSTART; x++) {
            #pragma ivdep
            #pragma vector always
            for ( y = YSTART; y < NY + YSTART; y++) {
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
