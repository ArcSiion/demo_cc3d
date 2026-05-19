#include "define.h"

void vectime_extra_array_unroll8(double* A, int NX, int NY, int NZ, int T) {
	if (back_scalar_256(NX, NY, NZ, T)) {
		naive_scalar(A, NX, NY, NZ, T);
		return;
	}
	
    double (* B)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;

    long int i, j, k;

    double tmp[4];
    int tt, t = 0, x, xx, y, yy, z, zz;

    vec v_y_plus_1, v_y_minus_1;
    vec v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7; 

    vec in, out;
	SET_COFF;

	double (* AV) [NY][NZ][VECLEN] = (double(*)[NY][NZ][VECLEN])alloc_extra_array(sizeof(double) * NY * NZ * VECLEN * 3);

	double (* BV0) [NZ][VECLEN] = (double (*) [NZ][VECLEN]) AV;
	double (* BV1) [NZ][VECLEN] = (double (*) [NZ][VECLEN]) (AV + 1);
    double (* BV2) [NZ][VECLEN] = (double (*) [NZ][VECLEN]) (AV + 2);

	double (* Btmp [3]) [NZ][VECLEN]  = {BV0, BV1, BV2};

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
				    transpose4x4(v_center_0, v_center_1, v_center_2, v_center_3, in, out, v_y_plus_1, v_y_minus_1);
                    vstore(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 0][0], v_center_0);
                    vstore(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 1][0], v_center_1);
                    vstore(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 2][0], v_center_2);
                    vstore(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 3][0], v_center_3);
                }
			}
		} 
        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VECLEN ; x ++){

            for ( y = YSTART; y < NY + YSTART; y++ ) {
                
                z = ZSTART;
                
                vset_3d(v_center_0, x, y, z-ZSLOPE);  // (0,0,-1)

                vload(v_center_1, BV1[y - YSTART][0][0]);

                for ( z = ZSTART; z <= NZ + ZSTART - VECLEN * 2; z += VECLEN * 2) {
                
                    vload(in, B[(t)%2][x+STRIDE*4][y][z]);   // the next x iter in vector

                    Compute_one(v_center_0, v_center_1, v_center_2, 1, 1);

                    Compute_one(v_center_1, v_center_2, v_center_3, 2, 2);

                    Compute_one(v_center_2, v_center_3, v_center_4, 3, 3);

                    Compute_one(v_center_3, v_center_4, v_center_5, 4, 4);

                    vstore(B[(t)%2][x][y][z], out);

                    vload(in, B[(t)%2][x+STRIDE*4][y][z + VECLEN]);   // the next x iter in vector

                    Compute_one(v_center_4, v_center_5, v_center_6, 5, 1);

                    Compute_one(v_center_5, v_center_6, v_center_7, 6, 2);

                    Compute_one(v_center_6, v_center_7, v_center_0, 7, 3);

                    Compute_last_one(v_center_7, v_center_0, v_center_1, 8, 4);
                    //-------------------------------------------------------------------------------------------

                    vstore(B[(t)%2][x][y][z + VECLEN], out);

                }
                if(z <= NZ + ZSTART - VECLEN) {
                
                    vload(in, B[(t)%2][x+STRIDE*4][y][z]);   // the next x iter in vector

                    Compute_one(v_center_0, v_center_1, v_center_2, 1, 1); 

                    Compute_one(v_center_1, v_center_2, v_center_3, 2, 2);

                    Compute_one(v_center_2, v_center_3, v_center_0, 3, 3);

 // there exists a trick of the v_center_1 loading,
 // the condition v_center_1 = (z > NZ + ZSTART - VECLEN - VECLEN * 2) is wrong
 // it should be  v_center_1 = (z > NZ + ZSTART - VECLEN - VECLEN)
 // but they are equal, they are both true.

                    Compute_last_one(v_center_3, v_center_0, v_center_1, 4, 4);

                    vstore(B[(t)%2][x][y][z], out);

                    z += VECLEN;
                }
                for ( ; z < NZ + ZSTART; z++) {

                    vset_3d(v_center_2 , x, y, z + 1);

                    Compute_1vector(v_center_1, v_center_0, v_center_2, \
                                    setv_3d(x, y - 1, z), \
                                    setv_3d(x, y + 1, z), \
                                    setv_3d(x - 1, y, z), \
                                    setv_3d(x + 1, y, z));
                    
                    vstore_set_3d(x, y, z, v_center_0);

                    v_center_0 = v_center_1;
                    v_center_1 = v_center_2;
                } 
            }
            Btmp[0] = BV0;
            BV0 = BV1;
            BV1 = BV2;
            BV2 = Btmp[0];     
        }
        Btmp [0] = BV0;
        Btmp [1] = BV1;
        Btmp [2] = BV2;
        ////when x  = NX + VECLEN - STRIDE * VECLEN,
        ////need extra transpose back to memory

        for(; x < NX + XSTART - STRIDE * VECLEN + 1 + 3; x++){           
			for ( y = YSTART ; y < NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VECLEN; z += VECLEN) {
                    vload(v_center_0, Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART)][(z - ZSTART) + 0][0]);
                    vload(v_center_1, Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART)][(z - ZSTART) + 1][0]);
                    vload(v_center_2, Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART)][(z - ZSTART) + 2][0]);
                    vload(v_center_3, Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART)][(z - ZSTART) + 3][0]);    
				    transpose4x4(v_center_0, v_center_1, v_center_2, v_center_3, in, out, v_y_plus_1, v_y_minus_1);
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
            #pragma ivdep
            #pragma vector always
            for ( z = 0; z < NZ + ZSTART * 2 ; z++) {
                B[1][XSTART - XSLOPE][y][z] = B[0][XSTART - XSLOPE][y][z];
            }
        }
	}
	//Extra points
	for ( ; t < T; t++){
		for (x = XSTART; x < NX + XSTART; x++) {
            for ( y = YSTART; y < NY + YSTART; y++) {
                #pragma ivdep
                #pragma vector always
                for( z = ZSTART; z < NZ + ZSTART; z++) {
                    Compute_scalar(B,t,x,y,z);
                }
            }
		}
	}	
    free_extra_array(AV);
}
