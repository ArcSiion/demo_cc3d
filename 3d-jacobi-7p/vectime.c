#include "define.h"

void vectime(double* A, int NX, int NY, int NZ, int T) {
	if (back_scalar_256(NX, NY, NZ, T)) {
		naive_scalar(A, NX, NY, NZ, T);
		return;
	}
	
    double (* B)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;

    long int i, j, k;

    double tmp[4];
    int tt, t =0, x, xx, y, yy, z, zz;

    vec v_x_plus_1, v_x_minus_1;
    vec v_y_plus_1, v_y_minus_1;
    vec v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5; 

    vec in, out;
	SET_COFF;

	for ( tt = 0; tt <= T - VECLEN; tt += VECLEN){	
		for( t = tt ; t < tt + VECLEN - 1 ; t++){		//head
			for ( x = XSTART; x < XSTART + STRIDE * (VECLEN - 1 - (t - tt)); x++) {//ASSERT VECLEN <= STRIDE + 1
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

        for(x = XSTART - XSLOPE; x <= XSTART +XSLOPE; x++){
			for ( y = YSTART; y < NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ; z += VECLEN ) {
				    vload(v_center_0, B[t%2    ][x + STRIDE *3][y][z]);
				    vload(v_center_1, B[(t+1)%2][x + STRIDE *2][y][z]);
				    vload(v_center_2, B[t%2    ][x + STRIDE *1][y][z]);
				    vload(v_center_3, B[(t+1)%2][x            ][y][z]);
				    transpose(v_center_0, v_center_1, v_center_2, v_center_3, in, out);
				    vstore( B[t%2    ][x + STRIDE *3][y][z], v_center_3);
				    vstore( B[(t+1)%2][x + STRIDE *2][y][z], v_center_2);
				    vstore( B[t%2    ][x + STRIDE *1][y][z], v_center_1);
				    vstore( B[(t+1)%2][x            ][y][z], v_center_0);
                }
			}
		} 
        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VECLEN ; x ++){
            for ( y = YSTART; y < NY + YSTART; y++ ) {
                z = ZSTART;
                vloadset(v_center_0, B[(t+1)%2][x         ][y][z-ZSLOPE], \
                                     B[(t)%2  ][x+STRIDE  ][y][z-ZSLOPE], \
                                     B[(t+1)%2][x+STRIDE*2][y][z-ZSLOPE], \
                                     B[(t)%2  ][x+STRIDE*3][y][z-ZSLOPE]);  // (0,0,-1)
                                    
                vload(v_center_1, B[(t+1)%2  ][x][y][z]);   //(0,0,0)
                for ( z = ZSTART; z <= NZ + ZSTART - VECLEN; z += VECLEN) {
                
                    vload(in, B[(t)%2][x+STRIDE*4][y][z]);   // the next x iter in vector

                                                                                        // (x, y,z)
                    vload(v_center_2 , B[(t)%2  ][x+STRIDE         ][y][z]);            // (0,0,+1)      

                    //&B[(t+1)%2][x  ][y-YSLOPE][z] (0,-1,0)
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y-YSLOPE][z], B[(t)%2][x+STRIDE][y-YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z], B[(t)%2][x+STRIDE*3][y-YSLOPE][z]))    \
                                  : _mm256_loadu_pd( &B[(t+1)%2][x  ][y-YSLOPE][z] );   // (0,-1,0)
                    //&B[(t+1)%2][x  ][y+YSLOPE][z] (0,+1,0)      
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y+YSLOPE][z], B[(t)%2][x+STRIDE][y+YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z], B[(t)%2][x+STRIDE*3][y+YSLOPE][z]))    \
                                  : _mm256_loadu_pd( &B[(t+1)%2][x  ][y+YSLOPE][z] );   // (0,+1,0)
                    vload(v_x_minus_1, B[(t+1)%2][x-XSLOPE         ][y][z]);            // (-1,0,0)
                    vload(v_x_plus_1 , B[(t+1)%2][x+XSLOPE         ][y][z]);            // (+1,0,0)

                    Compute_1vector(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //1st   store the newest value to the left vec

                    Input_Output_1(out, v_center_0, in);
                    vstore(B[(t+1)%2][x+STRIDE][y][z], v_center_0);   //for the compute of next x iteration
                    //-------------------------------------------------------------------------------------------

                    vload(v_center_3 , B[(t+1)%2][x+STRIDE*2       ][y][z]);                  // (0,0,+2)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y-YSLOPE][z+1], B[(t)%2][x+STRIDE][y-YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+1], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+1]))    \
                                  : _mm256_loadu_pd( &B[(t  )%2][x+STRIDE][y-YSLOPE][z] );    // (0,-1,+1)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y+YSLOPE][z+1], B[(t)%2][x+STRIDE][y+YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+1], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+1]))    \
                                  : _mm256_loadu_pd( &B[(t  )%2][x+STRIDE][y+YSLOPE][z] );    // (0,+1,+1)
                    vload(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE  ][y][z]);                  // (-1,0,+1)
                    vload(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE  ][y][z]);                  // (+1,0,+1)

                    Compute_1vector(v_center_2, v_center_1, v_center_3, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //2nd   store the newest value to the left vec

                    Input_Output_2(out, v_center_1, in);
                    vstore(B[(t)%2][x+STRIDE*2][y][z], v_center_1);  //for the compute of next x iteration
                    //-------------------------------------------------------------------------------------------

                    vload(v_center_4 , B[(t)%2  ][x+STRIDE*3       ][y][z]);                    // (0,0,+3)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y-YSLOPE][z+2], B[(t)%2][x+STRIDE][y-YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+2], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+2]))    \
                                  : _mm256_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z] );    // (0,-1,+2)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y+YSLOPE][z+2], B[(t)%2][x+STRIDE][y+YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+2], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+2]))    \
                                  : _mm256_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z] );    // (0,+1,+2)
                    vload(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*2][y][z]);                    // (-1,0,+2)
                    vload(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*2][y][z]);                    // (+1,0,+2)

                    Compute_1vector(v_center_3, v_center_2, v_center_4, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //3rd   store the newest value to the left vec

                    Input_Output_3(out, v_center_2, in);	
                    vstore(B[(t+1)%2  ][x+STRIDE*3][y][z], v_center_2);  //for the compute of next x iteration
                    //-------------------------------------------------------------------------------------------

                    v_center_5 = ( z > NZ + ZSTART - VECLEN - VECLEN ) ? \
                                 (_mm256_set_pd(B[(t+1)%2][x][y][z+VECLEN], B[(t)%2][x+STRIDE][y][z+VECLEN], \
                                                B[(t+1)%2][x+STRIDE*2][y][z+VECLEN], B[(t)%2][x+STRIDE*3][y][z+VECLEN])) \
                                : _mm256_loadu_pd(&B[(t+1)%2][x][y][z+VECLEN])  ;               // (0,0,+4) 
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y-YSLOPE][z+3], B[(t)%2][x+STRIDE][y-YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+3], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+3]))    \
                                  : _mm256_loadu_pd( &B[(t  )%2][x+STRIDE*3][y-YSLOPE][z] );    // (0,-1,+2)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm256_set_pd(B[(t+1)%2][x][y+YSLOPE][z+3], B[(t)%2][x+STRIDE][y+YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+3], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+3]))    \
                                  : _mm256_loadu_pd( &B[(t  )%2][x+STRIDE*3][y+YSLOPE][z] );    // (0,+1,+2)
                    vload(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*3][y][z]);                    // (-1,0,+3)
                    vload(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*3][y][z]);                    // (+1,0,+3)

                    Compute_1vector(v_center_4, v_center_3, v_center_5, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //4th   store the newest value to the left vec

                    Input_Output_4(out, v_center_3, in);	
                    vstore(B[(t)%2][x+STRIDE*4][y][z], v_center_3);  //for the compute of next x iteration
                    //-------------------------------------------------------------------------------------------

                    vstore(B[(t)%2][x][y][z], out);

                    v_center_0 = v_center_4;
                    v_center_1 = v_center_5;
                }
                for ( ; z < NZ + ZSTART; z++) {
                    vloadset(v_y_minus_1, B[(t+1)%2][x           ][y-1][z  ], \
                                          B[(t)%2  ][x+STRIDE    ][y-1][z  ], \
                                          B[(t+1)%2][x+STRIDE*2  ][y-1][z  ], \
                                          B[(t)%2  ][x+STRIDE*3  ][y-1][z  ]);
                    vloadset(v_y_plus_1 , B[(t+1)%2][x           ][y+1][z  ], \
                                          B[(t)%2  ][x+STRIDE    ][y+1][z  ], \
                                          B[(t+1)%2][x+STRIDE*2  ][y+1][z  ], \
                                          B[(t)%2  ][x+STRIDE*3  ][y+1][z  ]);
                    vloadset(v_x_minus_1, B[(t+1)%2][x-1         ][y  ][z  ], \
                                          B[(t)%2  ][x-1+STRIDE  ][y  ][z  ], \
                                          B[(t+1)%2][x-1+STRIDE*2][y  ][z  ], \
                                          B[(t)%2  ][x-1+STRIDE*3][y  ][z  ]);
                    vloadset(v_x_plus_1 , B[(t+1)%2][x+1         ][y  ][z  ], \
                                          B[(t)%2  ][x+1+STRIDE  ][y  ][z  ], \
                                          B[(t+1)%2][x+1+STRIDE*2][y  ][z  ], \
                                          B[(t)%2  ][x+1+STRIDE*3][y  ][z  ]);
                    vloadset(v_center_2 , B[(t+1)%2][x           ][y  ][z+1], \
                                          B[(t)%2  ][x+STRIDE    ][y  ][z+1], \
                                          B[(t+1)%2][x+STRIDE*2  ][y  ][z+1], \
                                          B[(t)%2  ][x+STRIDE*3  ][y  ][z+1]);

                    Compute_1vector(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    _mm256_storeu_pd(tmp, v_center_0);
                    B[(t+1)%2][x + STRIDE * 3][ y ][ z ] = tmp[0];
                    B[t%2    ][x + STRIDE * 2][ y ][ z ] = tmp[1];
                    B[(t+1)%2][x + STRIDE * 1][ y ][ z ] = tmp[2];
                    B[t%2    ][x             ][ y ][ z ] = tmp[3];

                    v_center_0 = v_center_1;
                    v_center_1 = v_center_2;
                } 
            }
        }
        ////when x  = NX + VECLEN - STRIDE * VECLEN,
        ////need extra transpose back to memory
        
        for(; x < NX + XSTART - STRIDE * VECLEN + 1 + 3; x++){           
			for ( y = YSTART ; y < NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VECLEN; z += VECLEN) {
				    vload(v_center_3, B[t%2    ][x - XSLOPE + STRIDE *3][y][z]);
				    vload(v_center_2, B[(t+1)%2][x - XSLOPE + STRIDE *2][y][z]);
				    vload(v_center_1, B[t%2    ][x - XSLOPE + STRIDE *1][y][z]);
				    vload(v_center_0, B[(t+1)%2][x - XSLOPE            ][y][z]);
				    transpose(v_center_0, v_center_1, v_center_2, v_center_3, in, out);
				    vstore( B[t%2    ][x - XSLOPE + STRIDE *3][y][z], v_center_0);
				    vstore( B[(t+1)%2][x - XSLOPE + STRIDE *2][y][z], v_center_1);
				    vstore( B[t%2    ][x - XSLOPE + STRIDE *1][y][z], v_center_2);
				    vstore( B[(t+1)%2][x - XSLOPE            ][y][z], v_center_3);
                }
			}
		}            

        ////tail
        xx = NX + XSTART - STRIDE * VECLEN + 1;
        for( t = tt ; t < tt + VECLEN ; t++){	
			for ( x = xx + STRIDE * (VECLEN - 1 - (t - tt)); x < NX + XSTART; x++) {
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
	//Extra points
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
