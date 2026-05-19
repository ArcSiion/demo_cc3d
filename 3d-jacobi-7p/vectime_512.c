#include "define_512.h"

void vectime_512(double* A, int NX, int NY, int NZ, int T) {
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
		for( t = tt ; t < tt + VVECLEN - 1 ; t++){		//head
			for ( x = XSTART; x < XSTART + STRIDE * (VVECLEN - 1 - (t - tt)); x++) {//ASSERT VVECLEN <= STRIDE + 1
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

        // 核心计算循环
        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VVECLEN ; x ++){
            for ( y = YSTART; y < NY + YSTART; y++ ) {
                z = ZSTART;
                
                vloadset_512(v_center_0, B[(t+1)%2][x         ][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE  ][y][z-ZSLOPE], \
                                         B[(t+1)%2][x+STRIDE*2][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE*3][y][z-ZSLOPE], \
                                         B[(t+1)%2][x+STRIDE*4][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE*5][y][z-ZSLOPE], \
                                         B[(t+1)%2][x+STRIDE*6][y][z-ZSLOPE], \
                                         B[(t)%2  ][x+STRIDE*7][y][z-ZSLOPE]);  // (0,0,-1)
                                    
                vload_512(v_center_1, B[(t+1)%2  ][x][y][z]);   //(0,0,0)
                
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                
                    vload_512(in, B[(t)%2][x+STRIDE*8][y][z]);   // the next x iter in vector

                    // =================================== 1st ===================================
                    vload_512(v_center_2 , B[(t)%2  ][x+STRIDE         ][y][z]);            // (0,0,+1)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z], B[(t)%2][x+STRIDE][y-YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z], B[(t)%2][x+STRIDE*3][y-YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z], B[(t)%2][x+STRIDE*5][y-YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z], B[(t)%2][x+STRIDE*7][y-YSLOPE][z]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x][y-YSLOPE][z] );   // (0,-1,0)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z], B[(t)%2][x+STRIDE][y+YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z], B[(t)%2][x+STRIDE*3][y+YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z], B[(t)%2][x+STRIDE*5][y+YSLOPE][z], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z], B[(t)%2][x+STRIDE*7][y+YSLOPE][z]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x  ][y+YSLOPE][z] );   // (0,+1,0)
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE         ][y][z]);            // (-1,0,0)
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE         ][y][z]);            // (+1,0,0)

                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //1st   store the newest value to the left vec
                    Input_Output_1_512(out, v_center_0, in);
                    vstore_512(B[(t+1)%2][x+STRIDE][y][z], v_center_0);   //for the compute of next x iteration
                    
                    // =================================== 2nd ===================================
                    vload_512(v_center_3 , B[(t+1)%2][x+STRIDE*2       ][y][z]);                  // (0,0,+2)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+1], B[(t)%2][x+STRIDE][y-YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+1], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+1], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+1], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+1]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y-YSLOPE][z] );    // (0,-1,+1)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+1], B[(t)%2][x+STRIDE][y+YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+1], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+1], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+1], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+1], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+1]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE][y+YSLOPE][z] );    // (0,+1,+1)
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE  ][y][z]);                  // (-1,0,+1)
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE  ][y][z]);                  // (+1,0,+1)

                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //2nd   
                    Input_Output_2_512(out, v_center_1, in);
                    vstore_512(B[(t)%2][x+STRIDE*2][y][z], v_center_1);  
                    
                    // =================================== 3rd ===================================
                    vload_512(v_center_4 , B[(t)%2  ][x+STRIDE*3       ][y][z]);                    // (0,0,+3)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+2], B[(t)%2][x+STRIDE][y-YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+2], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+2], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+2], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+2]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z] );    // (0,-1,+2)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+2], B[(t)%2][x+STRIDE][y+YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+2], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+2], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+2], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+2], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+2]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z] );    // (0,+1,+2)
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*2][y][z]);                    // (-1,0,+2)
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*2][y][z]);                    // (+1,0,+2)

                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //3rd   
                    Input_Output_3_512(out, v_center_2, in);	
                    vstore_512(B[(t+1)%2  ][x+STRIDE*3][y][z], v_center_2);  
                    
                    // =================================== 4th ===================================
                    vload_512(v_center_5 , B[(t+1)%2][x+STRIDE*4       ][y][z]);                  // (0,0,+4)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+3], B[(t)%2][x+STRIDE][y-YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+3], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+3], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+3], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+3]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y-YSLOPE][z] );    // (0,-1,+3)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+3], B[(t)%2][x+STRIDE][y+YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+3], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+3], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+3], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+3], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+3]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*3][y+YSLOPE][z] );    // (0,+1,+3)
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*3][y][z]);                  // (-1,0,+3)
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*3][y][z]);                  // (+1,0,+3)

                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //4th   
                    Input_Output_4_512(out, v_center_3, in);
                    vstore_512(B[(t)%2][x+STRIDE*4][y][z], v_center_3);  

                    // =================================== 5th ===================================
                    vload_512(v_center_6 , B[(t)%2  ][x+STRIDE*5       ][y][z]);                  // (0,0,+5)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+4], B[(t)%2][x+STRIDE][y-YSLOPE][z+4], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+4], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+4], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+4], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+4], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+4], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+4]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z] );    // (0,-1,+4)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+4], B[(t)%2][x+STRIDE][y+YSLOPE][z+4], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+4], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+4], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+4], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+4], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+4], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+4]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z] );    // (0,+1,+4)
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*4][y][z]);                  // (-1,0,+4)
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*4][y][z]);                  // (+1,0,+4)

                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //5th   
                    Input_Output_5_512(out, v_center_4, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*5][y][z], v_center_4);  

                    // =================================== 6th ===================================
                    vload_512(v_center_7 , B[(t+1)%2][x+STRIDE*6       ][y][z]);                  // (0,0,+6)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+5], B[(t)%2][x+STRIDE][y-YSLOPE][z+5], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+5], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+5], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+5], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+5], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+5], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+5]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y-YSLOPE][z] );    // (0,-1,+5)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+5], B[(t)%2][x+STRIDE][y+YSLOPE][z+5], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+5], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+5], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+5], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+5], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+5], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+5]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*5][y+YSLOPE][z] );    // (0,+1,+5)
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*5][y][z]);                  // (-1,0,+5)
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*5][y][z]);                  // (+1,0,+5)

                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //6th   
                    Input_Output_6_512(out, v_center_5, in);
                    vstore_512(B[(t)%2][x+STRIDE*6][y][z], v_center_5);  

                    // =================================== 7th ===================================
                    vload_512(v_center_8 , B[(t)%2  ][x+STRIDE*7       ][y][z]);                  // (0,0,+7)                                                                    
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+6], B[(t)%2][x+STRIDE][y-YSLOPE][z+6], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+6], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+6], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+6], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+6], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+6], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+6]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z] );    // (0,-1,+6)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+6], B[(t)%2][x+STRIDE][y+YSLOPE][z+6], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+6], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+6], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+6], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+6], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+6], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+6]))    \
                                  : _mm512_loadu_pd( &B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z] );    // (0,+1,+6)
                    vload_512(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*6][y][z]);                  // (-1,0,+6)
                    vload_512(v_x_plus_1 , B[(t+1)%2][x+XSLOPE+STRIDE*6][y][z]);                  // (+1,0,+6)

                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //7th   
                    Input_Output_7_512(out, v_center_6, in);
                    vstore_512(B[(t+1)%2][x+STRIDE*7][y][z], v_center_6);  

                    // =================================== 8th ===================================
                    v_center_9 = ( z > NZ + ZSTART - VVECLEN - VVECLEN ) ? \
                                 (_mm512_set_pd(B[(t+1)%2][x][y][z+VVECLEN], B[(t)%2][x+STRIDE][y][z+VVECLEN], \
                                                B[(t+1)%2][x+STRIDE*2][y][z+VVECLEN], B[(t)%2][x+STRIDE*3][y][z+VVECLEN], \
                                                B[(t+1)%2][x+STRIDE*4][y][z+VVECLEN], B[(t)%2][x+STRIDE*5][y][z+VVECLEN], \
                                                B[(t+1)%2][x+STRIDE*6][y][z+VVECLEN], B[(t)%2][x+STRIDE*7][y][z+VVECLEN])) \
                                : _mm512_loadu_pd(&B[(t+1)%2][x][y][z+VVECLEN])  ;               // (0,0,+8) 
                    v_y_minus_1 = ( y == YSTART ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y-YSLOPE][z+7], B[(t)%2][x+STRIDE][y-YSLOPE][z+7], \
                                                 B[(t+1)%2][x+STRIDE*2][y-YSLOPE][z+7], B[(t)%2][x+STRIDE*3][y-YSLOPE][z+7], \
                                                 B[(t+1)%2][x+STRIDE*4][y-YSLOPE][z+7], B[(t)%2][x+STRIDE*5][y-YSLOPE][z+7], \
                                                 B[(t+1)%2][x+STRIDE*6][y-YSLOPE][z+7], B[(t)%2][x+STRIDE*7][y-YSLOPE][z+7]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y-YSLOPE][z] );    // (0,-1,+7)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                  (_mm512_set_pd(B[(t+1)%2][x][y+YSLOPE][z+7], B[(t)%2][x+STRIDE][y+YSLOPE][z+7], \
                                                 B[(t+1)%2][x+STRIDE*2][y+YSLOPE][z+7], B[(t)%2][x+STRIDE*3][y+YSLOPE][z+7], \
                                                 B[(t+1)%2][x+STRIDE*4][y+YSLOPE][z+7], B[(t)%2][x+STRIDE*5][y+YSLOPE][z+7], \
                                                 B[(t+1)%2][x+STRIDE*6][y+YSLOPE][z+7], B[(t)%2][x+STRIDE*7][y+YSLOPE][z+7]))    \
                                  : _mm512_loadu_pd( &B[(t  )%2][x+STRIDE*7][y+YSLOPE][z] );    // (0,+1,+7)
                    vload_512(v_x_minus_1, B[(t  )%2][x-XSLOPE+STRIDE*7][y][z]);                    // (-1,0,+7)
                    vload_512(v_x_plus_1 , B[(t  )%2][x+XSLOPE+STRIDE*7][y][z]);                    // (+1,0,+7)

                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);   //8th   
                    Input_Output_8_512(out, v_center_7, in);	
                    vstore_512(B[(t)%2][x+STRIDE*8][y][z], v_center_7);  
                    
                    //-------------------------------------------------------------------------------------------

                    vstore_512(B[(t)%2][x][y][z], out);

                    // 滑动窗口推进 2 层 
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                }
                
                // 处理 Z 轴无法凑齐 8 个点结尾的不对齐区域
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

                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
                                    v_y_plus_1, v_x_minus_1, v_x_plus_1);
                    _mm512_storeu_pd(tmp, v_center_0);
                    B[(t+1)%2][x + STRIDE * 7][ y ][ z ] = tmp[0];
                    B[t%2    ][x + STRIDE * 6][ y ][ z ] = tmp[1];
                    B[(t+1)%2][x + STRIDE * 5][ y ][ z ] = tmp[2];
                    B[t%2    ][x + STRIDE * 4][ y ][ z ] = tmp[3];
                    B[(t+1)%2][x + STRIDE * 3][ y ][ z ] = tmp[4];
                    B[t%2    ][x + STRIDE * 2][ y ][ z ] = tmp[5];
                    B[(t+1)%2][x + STRIDE * 1][ y ][ z ] = tmp[6];
                    B[t%2    ][x             ][ y ][ z ] = tmp[7];

                    v_center_0 = v_center_1;
                    v_center_1 = v_center_2;
                } 
            }
        }
        
        // 转置，写回正确的内存槽位
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
