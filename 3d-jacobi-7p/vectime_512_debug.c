#include "define_512.h"

void vectime_512_debug(double* A, int NX, int NY, int NZ, int T) {
	
    double (* B)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;

    long int i, j, k;
    double tmp[VVECLEN];
    int tt, t =0, x, xx, y, yy, z, zz;

    vec_512 v_x_plus_1, v_x_minus_1;
    vec_512 v_y_plus_1, v_y_minus_1;
    vec_512 v_center_0, v_center_1, v_center_2, v_center_3, v_center_4;
    vec_512 v_center_5, v_center_6, v_center_7, v_center_8, v_center_9; 

    vec_512 in, out;
	SET_COFF_512;

	for ( tt = 0; tt <= T - VVECLEN; tt += VVECLEN){	
		for( t = tt ; t < tt + VVECLEN - 1 ; t++){		//head
			for ( x = XSTART; x < XSTART + STRIDE * (VVECLEN - 1 - (t - tt)); x++) {
                for ( y = YSTART; y < NY + YSTART; y++) {
                    #pragma ivdep
                    #pragma vector always
                    for ( z = ZSTART; z < NZ+ZSTART; z++) {
                        Compute_scalar(B, t, x, y, z);
                    }   
                }		
            }
		}
        t = tt;

        // 8x8 转置
        for(x = XSTART - XSLOPE; x <= XSTART +XSLOPE; x++){
			for ( y = YSTART; y < NY + YSTART ; y ++ ){
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
                    
				    vstore_512( B[t%2    ][x + STRIDE *7][y][z], v_center_7);
				    vstore_512( B[(t+1)%2][x + STRIDE *6][y][z], v_center_6);
				    vstore_512( B[t%2    ][x + STRIDE *5][y][z], v_center_5);
				    vstore_512( B[(t+1)%2][x + STRIDE *4][y][z], v_center_4);
				    vstore_512( B[t%2    ][x + STRIDE *3][y][z], v_center_3);
				    vstore_512( B[(t+1)%2][x + STRIDE *2][y][z], v_center_2);
				    vstore_512( B[t%2    ][x + STRIDE *1][y][z], v_center_1);
				    vstore_512( B[(t+1)%2][x            ][y][z], v_center_0);
                }
			}
		} 

        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VVECLEN ; x ++){
            y = YSTART;
            if (y < NY + YSTART) {
                z = ZSTART;
                vloadset_512(v_center_0, B[(t+1)%2][x         ][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE  ][y][z-ZSLOPE], \
                                         B[(t+1)%2][x+STRIDE*2][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE*3][y][z-ZSLOPE], \
                                         B[(t+1)%2][x+STRIDE*4][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE*5][y][z-ZSLOPE], \
                                         B[(t+1)%2][x+STRIDE*6][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE*7][y][z-ZSLOPE]);
                vload_512(v_center_1, B[(t+1)%2  ][x][y][z]);

                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]);

                    vload_512(v_center_2 , B[(t)%2  ][x+STRIDE][y][z]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z], B[(t)%2][x+STRIDE][y-YSLOPE][z], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z], B[(t)%2][x+STRIDE*3][y-YSLOPE][z], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z], B[(t)%2][x+STRIDE*5][y-YSLOPE][z], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z], B[(t)%2][x+STRIDE*7][y-YSLOPE][z]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x  ][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE][y][z]);
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE][y][z]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_1_512(out, v_center_0, in);
                    vstore_512(B[(t+1)%2][x+STRIDE][y][z], v_center_0);
                    
                    vload_512(v_center_3 , B[(t+1)%2][x+STRIDE*2][y][z]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+1], B[(t)%2][x+STRIDE][y-YSLOPE][z+1], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+1], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+1], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+1], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+1], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+1], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+1]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE][y][z]);
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE][y][z]);
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_2_512(out, v_center_1, in);
                    vstore_512(B[(t)%2][x+STRIDE*2][y][z], v_center_1);  
                    
                    vload_512(v_center_4 , B[(t)%2  ][x+STRIDE*3][y][z]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+2], B[(t)%2][x+STRIDE][y-YSLOPE][z+2], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+2], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+2], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+2], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+2], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+2], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+2]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*2][y][z]);
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*2][y][z]);
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_3_512(out, v_center_2, in);	
                    vstore_512(B[(t+1)%2  ][x+STRIDE*3][y][z], v_center_2);  
                    
                    vload_512(v_center_5 , B[(t+1)%2][x+STRIDE*4][y][z]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+3], B[(t)%2][x+STRIDE][y-YSLOPE][z+3], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+3], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+3], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+3], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+3], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+3], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+3]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*3][y][z]);
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*3][y][z]);
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_4_512(out, v_center_3, in);
                    vstore_512(B[(t)%2][x+STRIDE*4][y][z], v_center_3);  

                    vload_512(v_center_6 , B[(t)%2  ][x+STRIDE*5][y][z]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+4], B[(t)%2][x+STRIDE][y-YSLOPE][z+4], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+4], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+4], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+4], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+4], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+4], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+4]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*4][y][z]);
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*4][y][z]);
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_5_512(out, v_center_4, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*5][y][z], v_center_4);  

                    vload_512(v_center_7 , B[(t+1)%2][x+STRIDE*6][y][z]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+5], B[(t)%2][x+STRIDE][y-YSLOPE][z+5], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+5], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+5], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+5], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+5], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+5], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+5]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*5][y][z]);
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*5][y][z]);
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_6_512(out, v_center_5, in);
                    vstore_512(B[(t)%2][x+STRIDE*6][y][z], v_center_5);  

                    vload_512(v_center_8 , B[(t)%2  ][x+STRIDE*7][y][z]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+6], B[(t)%2][x+STRIDE][y-YSLOPE][z+6], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+6], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+6], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+6], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+6], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+6], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+6]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*6][y][z]);
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*6][y][z]);
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_7_512(out, v_center_6, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*7][y][z], v_center_6);  

                    v_center_9 = ( z > NZ + ZSTART - VVECLEN - VVECLEN ) ? (_mm512_set_pd(B[(t+1)%2][x][y][z+VVECLEN], B[(t)%2][x+STRIDE][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*2][y][z+VVECLEN], B[(t)%2][x+STRIDE*3][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*4][y][z+VVECLEN], B[(t)%2][x+STRIDE*5][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*6][y][z+VVECLEN], B[(t)%2][x+STRIDE*7][y][z+VVECLEN])) : _mm512_loadu_pd(&B[(t+1)%2][x][y][z+VVECLEN]);
                    v_y_minus_1 = _mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+7], B[(t)%2][x+STRIDE][y-YSLOPE][z+7], B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+7], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+7], B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+7], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+7], B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+7], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+7]);
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*7][y][z]);
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*7][y][z]);
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_8_512(out, v_center_7, in);	
                    vstore_512(B[(t)%2][x+STRIDE*8][y][z], v_center_7);  
                    
                    vstore_512(B[(t)%2][x][y][z], out);
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                }
                for ( ; z < NZ + ZSTART; z++) {
                    vloadset_512(v_y_minus_1, B[(t+1)%2][x           ][y-1][z  ], \
                                              B[(t)%2  ][x+STRIDE    ][y-1][z  ], \
                                              B[(t+1)%2][x+STRIDE*2  ][y-1][z  ], \
                                              B[(t)%2  ][x+STRIDE*3  ][y-1][z  ], \
                                              B[(t+1)%2][x+STRIDE*4  ][y-1][z  ], \
                                              B[(t)%2  ][x+STRIDE*5  ][y-1][z  ], \
                                              B[(t+1)%2][x+STRIDE*6  ][y-1][z  ], \
                                              B[(t)%2  ][x+STRIDE*7  ][y-1][z  ]);
                    vloadset_512(v_y_plus_1 , B[(t+1)%2][x           ][y+1][z  ], \
                                              B[(t)%2  ][x+STRIDE    ][y+1][z  ], \
                                              B[(t+1)%2][x+STRIDE*2  ][y+1][z  ], \
                                              B[(t)%2  ][x+STRIDE*3  ][y+1][z  ], \
                                              B[(t+1)%2][x+STRIDE*4  ][y+1][z  ], \
                                              B[(t)%2  ][x+STRIDE*5  ][y+1][z  ], \
                                              B[(t+1)%2][x+STRIDE*6  ][y+1][z  ], \
                                              B[(t)%2  ][x+STRIDE*7  ][y+1][z  ]);
                    vloadset_512(v_x_minus_1, B[(t+1)%2][x-1         ][y  ][z  ], \
                                              B[(t)%2  ][x-1+STRIDE  ][y  ][z  ], \
                                              B[(t+1)%2][x-1+STRIDE*2][y  ][z  ], \
                                              B[(t)%2  ][x-1+STRIDE*3][y  ][z  ], \
                                              B[(t+1)%2][x-1+STRIDE*4][y  ][z  ], \
                                              B[(t)%2  ][x-1+STRIDE*5][y  ][z  ], \
                                              B[(t+1)%2][x-1+STRIDE*6][y  ][z  ], \
                                              B[(t)%2  ][x-1+STRIDE*7][y  ][z  ]);
                    vloadset_512(v_x_plus_1 , B[(t+1)%2][x+1         ][y  ][z  ], \
                                              B[(t)%2  ][x+1+STRIDE  ][y  ][z  ], \
                                              B[(t+1)%2][x+1+STRIDE*2][y  ][z  ], \
                                              B[(t)%2  ][x+1+STRIDE*3][y  ][z  ], \
                                              B[(t+1)%2][x+1+STRIDE*4][y  ][z  ], \
                                              B[(t)%2  ][x+1+STRIDE*5][y  ][z  ], \
                                              B[(t+1)%2][x+1+STRIDE*6][y  ][z  ], \
                                              B[(t)%2  ][x+1+STRIDE*7][y  ][z  ]);
                    vloadset_512(v_center_2 , B[(t+1)%2][x           ][y  ][z+1], \
                                              B[(t)%2  ][x+STRIDE    ][y  ][z+1], \
                                              B[(t+1)%2][x+STRIDE*2  ][y  ][z+1], \
                                              B[(t)%2  ][x+STRIDE*3  ][y  ][z+1], \
                                              B[(t+1)%2][x+STRIDE*4  ][y  ][z+1], \
                                              B[(t)%2  ][x+STRIDE*5  ][y  ][z+1], \
                                              B[(t+1)%2][x+STRIDE*6  ][y  ][z+1], \
                                              B[(t)%2  ][x+STRIDE*7  ][y  ][z+1]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    _mm512_storeu_pd(tmp, v_center_0);
                    B[(t+1)%2][x+STRIDE*7][y][z]=tmp[0]; B[t%2][x+STRIDE*6][y][z]=tmp[1]; B[(t+1)%2][x+STRIDE*5][y][z]=tmp[2]; B[t%2][x+STRIDE*4][y][z]=tmp[3]; B[(t+1)%2][x+STRIDE*3][y][z]=tmp[4]; B[t%2][x+STRIDE*2][y][z]=tmp[5]; B[(t+1)%2][x+STRIDE*1][y][z]=tmp[6]; B[t%2][x][y][z]=tmp[7];
                    v_center_0 = v_center_1; v_center_1 = v_center_2;
                } 
            }
            for ( y = YSTART + 1; y < NY + YSTART - YSLOPE; y++ ) {
                z = ZSTART;
                vloadset_512(v_center_0, B[(t+1)%2][x][y][z-ZSLOPE], B[(t)%2][x+STRIDE][y][z-ZSLOPE], B[(t+1)%2][x+STRIDE*2][y][z-ZSLOPE], B[(t)%2][x+STRIDE*3][y][z-ZSLOPE], B[(t+1)%2][x+STRIDE*4][y][z-ZSLOPE], B[(t)%2][x+STRIDE*5][y][z-ZSLOPE], B[(t+1)%2][x+STRIDE*6][y][z-ZSLOPE], B[(t)%2][x+STRIDE*7][y][z-ZSLOPE]);
                vload_512(v_center_1, B[(t+1)%2][x][y][z]);   
                
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN - VVECLEN; z += VVECLEN) {
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]);

                    // 1st
                    vload_512(v_center_2 , B[(t)%2  ][x+STRIDE][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE][y][z]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_1_512(out, v_center_0, in);
                    vstore_512(B[(t+1)%2][x+STRIDE][y][z], v_center_0);
                    
                    // 2nd
                    vload_512(v_center_3 , B[(t+1)%2][x+STRIDE*2][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE][y][z]);
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_2_512(out, v_center_1, in);
                    vstore_512(B[(t)%2][x+STRIDE*2][y][z], v_center_1);  
                    
                    // 3rd
                    vload_512(v_center_4 , B[(t)%2  ][x+STRIDE*3][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*2][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*2][y][z]);
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_3_512(out, v_center_2, in);	
                    vstore_512(B[(t+1)%2  ][x+STRIDE*3][y][z], v_center_2);  
                    
                    // 4th
                    vload_512(v_center_5 , B[(t+1)%2][x+STRIDE*4][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*3][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*3][y][z]);
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_4_512(out, v_center_3, in);
                    vstore_512(B[(t)%2][x+STRIDE*4][y][z], v_center_3);  

                    // 5th
                    vload_512(v_center_6 , B[(t)%2  ][x+STRIDE*5][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*4][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*4][y][z]);
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_5_512(out, v_center_4, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*5][y][z], v_center_4);  

                    // 6th
                    vload_512(v_center_7 , B[(t+1)%2][x+STRIDE*6][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*5][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*5][y][z]);
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_6_512(out, v_center_5, in);
                    vstore_512(B[(t)%2][x+STRIDE*6][y][z], v_center_5);  

                    // 7th
                    vload_512(v_center_8 , B[(t)%2  ][x+STRIDE*7][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*6][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*6][y][z]);
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_7_512(out, v_center_6, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*7][y][z], v_center_6);  

                    // 8th (完全安全区)
                    v_center_9  = _mm512_loadu_pd(&B[(t+1)%2][x][y][z+VVECLEN]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*7][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*7][y][z]);
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_8_512(out, v_center_7, in);	
                    vstore_512(B[(t)%2][x+STRIDE*8][y][z], v_center_7);  
                    
                    vstore_512(B[(t)%2][x][y][z], out);
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                }
                
                if (z <= NZ + ZSTART - VVECLEN) {
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]);

                    vload_512(v_center_2 , B[(t)%2  ][x+STRIDE][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE][y][z]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_1_512(out, v_center_0, in); vstore_512(B[(t+1)%2][x+STRIDE][y][z], v_center_0);
                    
                    vload_512(v_center_3 , B[(t+1)%2][x+STRIDE*2][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE][y][z]);
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_2_512(out, v_center_1, in); vstore_512(B[(t)%2][x+STRIDE*2][y][z], v_center_1);  
                    
                    vload_512(v_center_4 , B[(t)%2  ][x+STRIDE*3][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*2][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*2][y][z]);
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_3_512(out, v_center_2, in); vstore_512(B[(t+1)%2  ][x+STRIDE*3][y][z], v_center_2);  
                    
                    vload_512(v_center_5 , B[(t+1)%2][x+STRIDE*4][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*3][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*3][y][z]);
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_4_512(out, v_center_3, in); vstore_512(B[(t)%2][x+STRIDE*4][y][z], v_center_3);  

                    vload_512(v_center_6 , B[(t)%2  ][x+STRIDE*5][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*4][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*4][y][z]);
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_5_512(out, v_center_4, in); vstore_512(B[(t+1)%2][x+STRIDE*5][y][z], v_center_4);  

                    vload_512(v_center_7 , B[(t+1)%2][x+STRIDE*6][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*5][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*5][y][z]);
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_6_512(out, v_center_5, in); vstore_512(B[(t)%2][x+STRIDE*6][y][z], v_center_5);  

                    vload_512(v_center_8 , B[(t)%2  ][x+STRIDE*7][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*6][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*6][y][z]);
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_7_512(out, v_center_6, in); vstore_512(B[(t+1)%2][x+STRIDE*7][y][z], v_center_6);  

                    v_center_9  = _mm512_set_pd(B[(t+1)%2][x][y][z+VVECLEN], B[(t)%2][x+STRIDE][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*2][y][z+VVECLEN], B[(t)%2][x+STRIDE*3][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*4][y][z+VVECLEN], B[(t)%2][x+STRIDE*5][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*6][y][z+VVECLEN], B[(t)%2][x+STRIDE*7][y][z+VVECLEN]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y-YSLOPE][z] ); v_y_plus_1  = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y+YSLOPE][z] );
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*7][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*7][y][z]);
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1); Input_Output_8_512(out, v_center_7, in);	
                    vstore_512(B[(t)%2][x+STRIDE*8][y][z], v_center_7);  
                    
                    vstore_512(B[(t)%2][x][y][z], out);
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                    z += VVECLEN;
                }
                
                for ( ; z < NZ + ZSTART; z++) {
                    vloadset_512(v_y_minus_1, B[(t+1)%2][x][y-1][z], B[(t)%2][x+STRIDE][y-1][z], B[(t+1)%2][x+STRIDE*2][y-1][z], B[(t)%2][x+STRIDE*3][y-1][z], B[(t+1)%2][x+STRIDE*4][y-1][z], B[(t)%2][x+STRIDE*5][y-1][z], B[(t+1)%2][x+STRIDE*6][y-1][z], B[(t)%2][x+STRIDE*7][y-1][z]);
                    vloadset_512(v_y_plus_1 , B[(t+1)%2][x][y+1][z], B[(t)%2][x+STRIDE][y+1][z], B[(t+1)%2][x+STRIDE*2][y+1][z], B[(t)%2][x+STRIDE*3][y+1][z], B[(t+1)%2][x+STRIDE*4][y+1][z], B[(t)%2][x+STRIDE*5][y+1][z], B[(t+1)%2][x+STRIDE*6][y+1][z], B[(t)%2][x+STRIDE*7][y+1][z]);
                    vloadset_512(v_x_minus_1, B[(t+1)%2][x-1][y][z], B[(t)%2][x-1+STRIDE][y][z], B[(t+1)%2][x-1+STRIDE*2][y][z], B[(t)%2][x-1+STRIDE*3][y][z], B[(t+1)%2][x-1+STRIDE*4][y][z], B[(t)%2][x-1+STRIDE*5][y][z], B[(t+1)%2][x-1+STRIDE*6][y][z], B[(t)%2][x-1+STRIDE*7][y][z]);
                    vloadset_512(v_x_plus_1 , B[(t+1)%2][x+1][y][z], B[(t)%2][x+1+STRIDE][y][z], B[(t+1)%2][x+1+STRIDE*2][y][z], B[(t)%2][x+1+STRIDE*3][y][z], B[(t+1)%2][x+1+STRIDE*4][y][z], B[(t)%2][x+1+STRIDE*5][y][z], B[(t+1)%2][x+1+STRIDE*6][y][z], B[(t)%2][x+1+STRIDE*7][y][z]);
                    vloadset_512(v_center_2 , B[(t+1)%2][x][y][z+1], B[(t)%2][x+STRIDE][y][z+1], B[(t+1)%2][x+STRIDE*2][y][z+1], B[(t)%2][x+STRIDE*3][y][z+1], B[(t+1)%2][x+STRIDE*4][y][z+1], B[(t)%2][x+STRIDE*5][y][z+1], B[(t+1)%2][x+STRIDE*6][y][z+1], B[(t)%2][x+STRIDE*7][y][z+1]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    _mm512_storeu_pd(tmp, v_center_0);
                    B[(t+1)%2][x+STRIDE*7][y][z]=tmp[0]; B[t%2][x+STRIDE*6][y][z]=tmp[1]; B[(t+1)%2][x+STRIDE*5][y][z]=tmp[2]; B[t%2][x+STRIDE*4][y][z]=tmp[3]; B[(t+1)%2][x+STRIDE*3][y][z]=tmp[4]; B[t%2][x+STRIDE*2][y][z]=tmp[5]; B[(t+1)%2][x+STRIDE*1][y][z]=tmp[6]; B[t%2][x][y][z]=tmp[7];
                    v_center_0 = v_center_1; v_center_1 = v_center_2;
                } 
            }

            int y_bottom = (YSTART + 1 > NY + YSTART - YSLOPE) ? YSTART + 1 : NY + YSTART - YSLOPE;
            for ( y = y_bottom; y < NY + YSTART; y++ ) {
                z = ZSTART;
                vloadset_512(v_center_0, B[(t+1)%2][x         ][y][z-ZSLOPE], B[(t)%2  ][x+STRIDE  ][y][z-ZSLOPE], B[(t+1)%2][x+STRIDE*2][y][z-ZSLOPE], B[(t)%2  ][x+STRIDE*3][y][z-ZSLOPE], B[(t+1)%2][x+STRIDE*4][y][z-ZSLOPE], B[(t)%2  ][x+STRIDE*5][y][z-ZSLOPE], B[(t+1)%2][x+STRIDE*6][y][z-ZSLOPE], B[(t)%2  ][x+STRIDE*7][y][z-ZSLOPE]);
                vload_512(v_center_1, B[(t+1)%2  ][x][y][z]);

                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]);

                    vload_512(v_center_2 , B[(t)%2  ][x+STRIDE][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z], B[(t)%2][x+STRIDE][y+YSLOPE][z], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z], B[(t)%2][x+STRIDE*3][y+YSLOPE][z], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z], B[(t)%2][x+STRIDE*5][y+YSLOPE][z], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z], B[(t)%2][x+STRIDE*7][y+YSLOPE][z]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE][y][z]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_1_512(out, v_center_0, in);
                    vstore_512(B[(t+1)%2][x+STRIDE][y][z], v_center_0);
                    
                    vload_512(v_center_3 , B[(t+1)%2][x+STRIDE*2][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+1], B[(t)%2][x+STRIDE][y+YSLOPE][z+1], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+1], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+1], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+1], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+1], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+1], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+1]);
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE][y][z]);
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_2_512(out, v_center_1, in);
                    vstore_512(B[(t)%2][x+STRIDE*2][y][z], v_center_1);  
                    
                    vload_512(v_center_4 , B[(t)%2  ][x+STRIDE*3][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+2], B[(t)%2][x+STRIDE][y+YSLOPE][z+2], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+2], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+2], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+2], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+2], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+2], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+2]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*2][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*2][y][z]);
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_3_512(out, v_center_2, in);	
                    vstore_512(B[(t+1)%2  ][x+STRIDE*3][y][z], v_center_2);  
                    
                    vload_512(v_center_5 , B[(t+1)%2][x+STRIDE*4][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+3], B[(t)%2][x+STRIDE][y+YSLOPE][z+3], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+3], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+3], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+3], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+3], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+3], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+3]);
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*3][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*3][y][z]);
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_4_512(out, v_center_3, in);
                    vstore_512(B[(t)%2][x+STRIDE*4][y][z], v_center_3);  

                    vload_512(v_center_6 , B[(t)%2  ][x+STRIDE*5][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+4], B[(t)%2][x+STRIDE][y+YSLOPE][z+4], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+4], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+4], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+4], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+4], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+4], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+4]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*4][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*4][y][z]);
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_5_512(out, v_center_4, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*5][y][z], v_center_4);  

                    vload_512(v_center_7 , B[(t+1)%2][x+STRIDE*6][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+5], B[(t)%2][x+STRIDE][y+YSLOPE][z+5], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+5], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+5], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+5], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+5], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+5], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+5]);
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*5][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*5][y][z]);
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_6_512(out, v_center_5, in);
                    vstore_512(B[(t)%2][x+STRIDE*6][y][z], v_center_5);  

                    vload_512(v_center_8 , B[(t)%2  ][x+STRIDE*7][y][z]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+6], B[(t)%2][x+STRIDE][y+YSLOPE][z+6], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+6], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+6], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+6], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+6], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+6], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+6]);
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*6][y][z]); vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*6][y][z]);
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_7_512(out, v_center_6, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*7][y][z], v_center_6);  

                    v_center_9 = ( z > NZ + ZSTART - VVECLEN - VVECLEN ) ? (_mm512_set_pd(B[(t+1)%2][x][y][z+VVECLEN], B[(t)%2][x+STRIDE][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*2][y][z+VVECLEN], B[(t)%2][x+STRIDE*3][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*4][y][z+VVECLEN], B[(t)%2][x+STRIDE*5][y][z+VVECLEN], B[(t+1)%2][x+STRIDE*6][y][z+VVECLEN], B[(t)%2][x+STRIDE*7][y][z+VVECLEN])) : _mm512_loadu_pd(&B[(t+1)%2][x][y][z+VVECLEN]);
                    v_y_minus_1 = _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y-YSLOPE][z] );
                    v_y_plus_1  = _mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+7], B[(t)%2][x+STRIDE][y+YSLOPE][z+7], B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+7], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+7], B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+7], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+7], B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+7], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+7]);
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*7][y][z]); vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*7][y][z]);
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    Input_Output_8_512(out, v_center_7, in);	
                    vstore_512(B[(t)%2][x+STRIDE*8][y][z], v_center_7);  
                    
                    vstore_512(B[(t)%2][x][y][z], out);
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                }
                for ( ; z < NZ + ZSTART; z++) {
                    vloadset_512(v_y_minus_1, B[(t+1)%2][x][y-1][z], B[(t)%2][x+STRIDE][y-1][z], B[(t+1)%2][x+STRIDE*2][y-1][z], B[(t)%2][x+STRIDE*3][y-1][z], B[(t+1)%2][x+STRIDE*4][y-1][z], B[(t)%2][x+STRIDE*5][y-1][z], B[(t+1)%2][x+STRIDE*6][y-1][z], B[(t)%2][x+STRIDE*7][y-1][z]);
                    vloadset_512(v_y_plus_1 , B[(t+1)%2][x][y+1][z], B[(t)%2][x+STRIDE][y+1][z], B[(t+1)%2][x+STRIDE*2][y+1][z], B[(t)%2][x+STRIDE*3][y+1][z], B[(t+1)%2][x+STRIDE*4][y+1][z], B[(t)%2][x+STRIDE*5][y+1][z], B[(t+1)%2][x+STRIDE*6][y+1][z], B[(t)%2][x+STRIDE*7][y+1][z]);
                    vloadset_512(v_x_minus_1, B[(t+1)%2][x-1][y][z], B[(t)%2][x-1+STRIDE][y][z], B[(t+1)%2][x-1+STRIDE*2][y][z], B[(t)%2][x-1+STRIDE*3][y][z], B[(t+1)%2][x-1+STRIDE*4][y][z], B[(t)%2][x-1+STRIDE*5][y][z], B[(t+1)%2][x-1+STRIDE*6][y][z], B[(t)%2][x-1+STRIDE*7][y][z]);
                    vloadset_512(v_x_plus_1 , B[(t+1)%2][x+1][y][z], B[(t)%2][x+1+STRIDE][y][z], B[(t+1)%2][x+1+STRIDE*2][y][z], B[(t)%2][x+1+STRIDE*3][y][z], B[(t+1)%2][x+1+STRIDE*4][y][z], B[(t)%2][x+1+STRIDE*5][y][z], B[(t+1)%2][x+1+STRIDE*6][y][z], B[(t)%2][x+1+STRIDE*7][y][z]);
                    vloadset_512(v_center_2 , B[(t+1)%2][x][y][z+1], B[(t)%2][x+STRIDE][y][z+1], B[(t+1)%2][x+STRIDE*2][y][z+1], B[(t)%2][x+STRIDE*3][y][z+1], B[(t+1)%2][x+STRIDE*4][y][z+1], B[(t)%2][x+STRIDE*5][y][z+1], B[(t+1)%2][x+STRIDE*6][y][z+1], B[(t)%2][x+STRIDE*7][y][z+1]);
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    _mm512_storeu_pd(tmp, v_center_0);
                    B[(t+1)%2][x+STRIDE*7][y][z]=tmp[0]; B[t%2][x+STRIDE*6][y][z]=tmp[1]; B[(t+1)%2][x+STRIDE*5][y][z]=tmp[2]; B[t%2][x+STRIDE*4][y][z]=tmp[3]; B[(t+1)%2][x+STRIDE*3][y][z]=tmp[4]; B[t%2][x+STRIDE*2][y][z]=tmp[5]; B[(t+1)%2][x+STRIDE*1][y][z]=tmp[6]; B[t%2][x][y][z]=tmp[7];
                    v_center_0 = v_center_1; v_center_1 = v_center_2;
                } 
            }
        }
        
        // 扫尾转置，将倾斜的寄存器波前写回正确的内存槽位
        for(; x < NX + XSTART - STRIDE * VVECLEN + 1 + 3; x++){           
			for ( y = YSTART ; y < NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
				    vload_512(v_center_7, B[t%2    ][x - XSLOPE + STRIDE *7][y][z]);
				    vload_512(v_center_6, B[(t+1)%2][x - XSLOPE + STRIDE *6][y][z]);
				    vload_512(v_center_5, B[t%2    ][x - XSLOPE + STRIDE *5][y][z]);
				    vload_512(v_center_4, B[(t+1)%2][x - XSLOPE + STRIDE *4][y][z]);
				    vload_512(v_center_3, B[t%2    ][x - XSLOPE + STRIDE *3][y][z]);
				    vload_512(v_center_2, B[(t+1)%2][x - XSLOPE + STRIDE *2][y][z]);
				    vload_512(v_center_1, B[t%2    ][x - XSLOPE + STRIDE *1][y][z]);
				    vload_512(v_center_0, B[(t+1)%2][x - XSLOPE            ][y][z]);
                    
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

        // Tail 残余区域标量处理
        xx = NX + XSTART - STRIDE * VVECLEN + 1;
        for( t = tt ; t < tt + VVECLEN ; t++){	
			for ( x = xx + STRIDE * (VVECLEN - 1 - (t - tt)); x < NX + XSTART; x++) {
                for ( y = YSTART; y < NY + YSTART; y++) {
                    #pragma ivdep
                    #pragma vector always
                    for ( z = ZSTART; z < NZ + ZSTART; z++) {
                        Compute_scalar(B,t,x,y,z);
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
    
	// 尾部无法填满 VVECLEN 时间块 Extra points
	for ( ; t < T; t++){
		for (x = XSTART; x < NX + XSTART; x++) {
            #pragma ivdep
            #pragma vector always
            for ( y = YSTART; y < NY + YSTART; y++) {
                for( z = ZSTART; z < NZ + ZSTART; z++) {
                    Compute_scalar(B,t,x,y,z);
                }
            }
		}
	}	
}