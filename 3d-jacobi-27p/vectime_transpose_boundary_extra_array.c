#include "define.h"

void vectime_transpose_boundary_extra_array(double* A, int NX, int NY, int NZ, int T) {
	if (back_scalar_256(NX, NY, NZ, T)) {
		naive_scalar(A, NX, NY, NZ, T);
		return;
	}
	
    double (* B)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;

    long int i, j, k;

    double tmp[VECLEN];
    int tt, t =0, x, xx, y, yy, z, zz;
    
    vec v_center_0, v_center_1, v_center_2, v_center_3;
    vec v_all_d_1_0, v_all_d_1_1, v_all_d_1_2, v_all_d_1_3;
    vec v_all_d_2_0, v_all_d_2_1, v_all_d_2_2, v_all_d_2_3;

    vec in, out;
	SET_COFF;

	double (* AV) [NY + 2][NZ][VECLEN] = (double(*)[NY+2][NZ][VECLEN])alloc_extra_array(sizeof(double) * (NY + 2) * NZ * VECLEN * 4);

	double (* BV0) [NZ][VECLEN] = (double (*) [NZ][VECLEN]) AV;
	double (* BV1) [NZ][VECLEN] = (double (*) [NZ][VECLEN]) (AV + 1);
    double (* BV2) [NZ][VECLEN] = (double (*) [NZ][VECLEN]) (AV + 2);
    double (* BV3) [NZ][VECLEN] = (double (*) [NZ][VECLEN]) (AV + 3);

	double (* Btmp [4]) [NZ][VECLEN]  = {BV0, BV1, BV2, BV3};

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
			for ( y = YSTART - YSLOPE; y <= NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ; z += VECLEN ) {
				    vload(v_center_0, B[t%2    ][x + STRIDE *3][y][z]);
				    vload(v_center_1, B[(t+1)%2][x + STRIDE *2][y][z]);
				    vload(v_center_2, B[t%2    ][x + STRIDE *1][y][z]);
				    vload(v_center_3, B[(t+1)%2][x            ][y][z]);
				    transpose(v_center_0, v_center_1, v_center_2, v_center_3, in, out);
                    _mm256_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 0][0], v_center_0);
                    _mm256_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 1][0], v_center_1);
                    _mm256_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 2][0], v_center_2);
                    _mm256_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 3][0], v_center_3);
                }
			}
		}

        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VECLEN ; x ++){
 
            for ( y = YSTART; y < NY + YSTART; y++ ) {

                z = ZSTART;

                v_center_0 =  load_v2(x, y, z, 0, 0, -1);
                v_all_d_1_0 = Add_4_vectors(    load_v2(x, y, z, 1, 0, -1),\
                                                load_v2(x, y, z, 0, 1, -1),\
                                                load_v2(x, y, z, -1, 0, -1),\
                                                load_v2(x, y, z, 0, -1, -1));
                v_all_d_2_0 = Add_4_vectors(    load_v2(x, y, z, 1, 1, -1),\
                                                load_v2(x, y, z, -1, 1, -1),\
                                                load_v2(x, y, z, -1, -1, -1),\
                                                load_v2(x, y, z, 1, -1, -1));

                v_center_1 = _mm256_loadu_pd( & BV1[y - YSTART + YSLOPE][0][0]);
                v_all_d_1_1 = Add_4_vectors(    load_x_m(x, y, z, 0, 0),\
                                                load_x_p(x, y, z, 0, 0),\
                                                load_x_c(x, y, z, 1, 0),\
                                                load_x_c(x, y, z, -1, 0));
                v_all_d_2_1 = Add_4_vectors(    load_x_m(x, y, z, 1, 0),\
                                                load_x_m(x, y, z, -1, 0),\
                                                load_x_p(x, y, z, 1, 0),\
                                                load_x_p(x, y, z, -1, 0));


                for ( z = ZSTART; z <= NZ + ZSTART - VECLEN; z += VECLEN) {
                
                    vload(in, B[(t)%2][x+STRIDE*VECLEN][y][z]);   // the next x iter in vector

                    v_center_2 = load_x_c(x, y, z, 0, 1); //_mm256_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 1][0]);
                    v_all_d_1_2 = Add_4_vectors(    load_x_m(x, y, z, 0, 1),\
                                                    load_x_p(x, y, z, 0, 1),\
                                                    load_x_c(x, y, z, 1, 1),\
                                                    load_x_c(x, y, z, -1, 1));
                    v_all_d_2_2 = Add_4_vectors(    load_x_m(x, y, z, 1, 1),\
                                                    load_x_m(x, y, z, -1, 1),\
                                                    load_x_p(x, y, z, 1, 1),\
                                                    load_x_p(x, y, z, -1, 1));
                    Compute_1vector(v_center_0, \
                                    v_center_1, \
                                    v_center_2, \
                                    v_all_d_1_0, \
                                    v_all_d_1_1, \
                                    v_all_d_1_2, \
                                    v_all_d_2_0, \
                                    v_all_d_2_1, \
                                    v_all_d_2_2);
                    Input_Output_1(out, v_center_0, in);
                    vstore(BV3[y - YSTART + YSLOPE][z - ZSTART][0], v_center_0);
                    //-------------------------------------------------------------------------------------------
                                                           
                    v_center_3 = _mm256_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 2][0]);
                    v_all_d_1_3 = Add_4_vectors(    load_x_m(x, y, z, 0, 2),\
                                                    load_x_p(x, y, z, 0, 2),\
                                                    load_x_c(x, y, z, 1, 2),\
                                                    load_x_c(x, y, z, -1, 2));
                    v_all_d_2_3 = Add_4_vectors(    load_x_m(x, y, z, 1, 2),\
                                                    load_x_m(x, y, z, -1, 2),\
                                                    load_x_p(x, y, z, 1, 2),\
                                                    load_x_p(x, y, z, -1, 2));
                    Compute_1vector(v_center_1, \
                                    v_center_2, \
                                    v_center_3, \
                                    v_all_d_1_1, \
                                    v_all_d_1_2, \
                                    v_all_d_1_3, \
                                    v_all_d_2_1, \
                                    v_all_d_2_2, \
                                    v_all_d_2_3);
                    Input_Output_2(out, v_center_1, in);
                    vstore(BV3[y - YSTART + YSLOPE][z - ZSTART + 1][0], v_center_1);
                    //-------------------------------------------------------------------------------------------
                                                                
                    v_center_0 = _mm256_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 3][0]);
                    v_all_d_1_0 = Add_4_vectors(    load_x_m(x, y, z, 0, 3),\
                                                    load_x_p(x, y, z, 0, 3),\
                                                    load_x_c(x, y, z, 1, 3),\
                                                    load_x_c(x, y, z, -1, 3));
                    v_all_d_2_0 = Add_4_vectors(    load_x_m(x, y, z, 1, 3),\
                                                    load_x_m(x, y, z, -1, 3),\
                                                    load_x_p(x, y, z, 1, 3),\
                                                    load_x_p(x, y, z, -1, 3));
                    Compute_1vector(v_center_2, \
                                    v_center_3, \
                                    v_center_0, \
                                    v_all_d_1_2, \
                                    v_all_d_1_3, \
                                    v_all_d_1_0, \
                                    v_all_d_2_2, \
                                    v_all_d_2_3, \
                                    v_all_d_2_0);
                    Input_Output_3(out, v_center_2, in);	
                    vstore(BV3[y - YSTART + YSLOPE][z - ZSTART + 2][0], v_center_2);
                    //-------------------------------------------------------------------------------------------
                    
                    if ( z > NZ + ZSTART - VECLEN - VECLEN ){
                        v_center_1 =  load_v2(x, y, z, 0, 0, VECLEN);
                        v_all_d_1_1 = Add_4_vectors(    load_v2(x, y, z, 1, 0, VECLEN),\
                                                        load_v2(x, y, z, 0, 1, VECLEN),\
                                                        load_v2(x, y, z, -1, 0, VECLEN),\
                                                        load_v2(x, y, z, 0, -1, VECLEN));
                        v_all_d_2_1 = Add_4_vectors(    load_v2(x, y, z, 1, 1, VECLEN),\
                                                        load_v2(x, y, z, -1, 1, VECLEN),\
                                                        load_v2(x, y, z, -1, -1, VECLEN),\
                                                        load_v2(x, y, z, 1, -1, VECLEN));
                    } else{
                        v_center_1 = _mm256_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 4][0]);
                        v_all_d_1_1 = Add_4_vectors(    load_x_m(x, y, z, 0, 4),\
                                                        load_x_p(x, y, z, 0, 4),\
                                                        load_x_c(x, y, z, 1, 4),\
                                                        load_x_c(x, y, z, -1, 4));
                        v_all_d_2_1 = Add_4_vectors(    load_x_m(x, y, z, 1, 4),\
                                                        load_x_m(x, y, z, -1, 4),\
                                                        load_x_p(x, y, z, 1, 4),\
                                                        load_x_p(x, y, z, -1, 4));
                    }

                    Compute_1vector(v_center_3, \
                                    v_center_0, \
                                    v_center_1, \
                                    v_all_d_1_3, \
                                    v_all_d_1_0, \
                                    v_all_d_1_1, \
                                    v_all_d_2_3, \
                                    v_all_d_2_0, \
                                    v_all_d_2_1);
                    Input_Output_4(out, v_center_3, in);	
                    vstore(BV3[y - YSTART + YSLOPE][z - ZSTART + 3][0], v_center_3);
                    //-------------------------------------------------------------------------------------------

                    vstore(B[(t)%2][x][y][z], out);

                }
                // for ( ; z < NZ + ZSTART; z++) {
                //     vloadset(v_y_minus_1, B[(t+1)%2][x           ][y-1][z  ], \
                //                           B[(t)%2  ][x+STRIDE    ][y-1][z  ], \
                //                           B[(t+1)%2][x+STRIDE*2  ][y-1][z  ], \
                //                           B[(t)%2  ][x+STRIDE*3  ][y-1][z  ]);
                //     vloadset(v_y_plus_1 , B[(t+1)%2][x           ][y+1][z  ], \
                //                           B[(t)%2  ][x+STRIDE    ][y+1][z  ], \
                //                           B[(t+1)%2][x+STRIDE*2  ][y+1][z  ], \
                //                           B[(t)%2  ][x+STRIDE*3  ][y+1][z  ]);
                //     vloadset(v_x_minus_1, B[(t+1)%2][x-1         ][y  ][z  ], \
                //                           B[(t)%2  ][x-1+STRIDE  ][y  ][z  ], \
                //                           B[(t+1)%2][x-1+STRIDE*2][y  ][z  ], \
                //                           B[(t)%2  ][x-1+STRIDE*3][y  ][z  ]);
                //     vloadset(v_x_plus_1 , B[(t+1)%2][x+1         ][y  ][z  ], \
                //                           B[(t)%2  ][x+1+STRIDE  ][y  ][z  ], \
                //                           B[(t+1)%2][x+1+STRIDE*2][y  ][z  ], \
                //                           B[(t)%2  ][x+1+STRIDE*3][y  ][z  ]);
                //     vloadset(v_center_2 , B[(t+1)%2][x           ][y  ][z+1], \
                //                           B[(t)%2  ][x+STRIDE    ][y  ][z+1], \
                //                           B[(t+1)%2][x+STRIDE*2  ][y  ][z+1], \
                //                           B[(t)%2  ][x+STRIDE*3  ][y  ][z+1]);

                //     Compute_1vector(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
                //                     v_y_plus_1, v_x_minus_1, v_x_plus_1);
                //     _mm256_storeu_pd(tmp, v_center_0);
                //     B[(t+1)%2][x + STRIDE * 3][ y ][ z ] = tmp[0];
                //     B[t%2    ][x + STRIDE * 2][ y ][ z ] = tmp[1];
                //     B[(t+1)%2][x + STRIDE * 1][ y ][ z ] = tmp[2];
                //     B[t%2    ][x             ][ y ][ z ] = tmp[3];

                //     v_center_0 = v_center_1;
                //     v_center_1 = v_center_2;
                // }                                
                for ( ; z < NZ + ZSTART; z++) {

                    v_center_2 = setv_3d(x, y, z + 1);
                    v_all_d_1_2 = Add_4_vectors(
                        setv_3d(x + 1, y    , z + 1),
                        setv_3d(x    , y + 1, z + 1),
                        setv_3d(x - 1, y    , z + 1),
                        setv_3d(x    , y - 1, z + 1)
                    );
                    v_all_d_2_2 = Add_4_vectors(
                        setv_3d(x + 1, y + 1, z + 1),
                        setv_3d(x - 1, y + 1, z + 1),
                        setv_3d(x - 1, y - 1, z + 1),
                        setv_3d(x + 1, y - 1, z + 1)
                    );

                    Compute_1vector(
                        v_center_0,
                        v_center_1,
                        v_center_2,
                        v_all_d_1_0,
                        v_all_d_1_1,
                        v_all_d_1_2,
                        v_all_d_2_0,
                        v_all_d_2_1,
                        v_all_d_2_2
                    );

                    vstore_set_3d(x, y, z, v_center_0);

                    v_center_0  = v_center_1;
                    v_center_1  = v_center_2;
                    v_all_d_1_0 = v_all_d_1_1;
                    v_all_d_1_1 = v_all_d_1_2;
                    v_all_d_2_0 = v_all_d_2_1;
                    v_all_d_2_1 = v_all_d_2_2;
                }
            }

            y = YSTART - YSLOPE;
			for ( z = ZSTART ; z <= NZ + ZSTART - VECLEN; z += VECLEN){	
                vload(v_center_0, B[t%2    ][x + XSLOPE + 1 + STRIDE *3][YSTART - YSLOPE][z]);
				vload(v_center_1, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][YSTART - YSLOPE][z]);
				vload(v_center_2, B[t%2    ][x + XSLOPE + 1 + STRIDE *1][YSTART - YSLOPE][z]);
				vload(v_center_3, B[(t+1)%2][x + XSLOPE + 1            ][YSTART - YSLOPE][z]);

                vload(v_all_d_1_0 , B[t%2    ][x + XSLOPE + 1 + STRIDE *3][YSTART + NY][z]); //y = NY + YSTART
				vload(v_all_d_1_1, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][YSTART + NY][z]);
				vload(v_all_d_1_2 , B[t%2    ][x + XSLOPE + 1 + STRIDE *1][YSTART + NY][z]);
				vload(v_all_d_1_3, B[(t+1)%2][x + XSLOPE + 1            ][YSTART + NY][z]);

				transpose(v_center_0, v_center_1, v_center_2, v_center_3, in, out);
				transpose(v_all_d_1_0, v_all_d_1_1, v_all_d_1_2, v_all_d_1_3, in, out);   //y = NY + YSTART

                _mm256_storeu_pd(&BV3[0][(z - ZSTART) + 3][0], v_center_3);
                _mm256_storeu_pd(&BV3[0][(z - ZSTART) + 2][0], v_center_2);
                _mm256_storeu_pd(&BV3[0][(z - ZSTART) + 1][0], v_center_1);
                _mm256_storeu_pd(&BV3[0][(z - ZSTART) + 0][0], v_center_0);

                _mm256_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 0][0], v_all_d_1_0);
                _mm256_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 1][0], v_all_d_1_1);
                _mm256_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 2][0], v_all_d_1_2);
                _mm256_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 3][0], v_all_d_1_3);
			}
            Btmp[0] = BV0;
            BV0 = BV1;
            BV1 = BV2;
            BV2 = BV3;
            BV3 = Btmp[0];   
        }
        Btmp [0] = BV0;
        Btmp [1] = BV1;
        Btmp [2] = BV2;
        Btmp [3] = BV3;

        for(; x < NX + XSTART - STRIDE * VECLEN + 1 + 3; x++){           
			for ( y = YSTART - YSLOPE ; y <= NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VECLEN; z += VECLEN) {
                    v_center_0 = _mm256_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 0][0]);
                    v_center_1 = _mm256_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 1][0]);
                    v_center_2 = _mm256_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 2][0]);
                    v_center_3 = _mm256_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 3][0]);  
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
    free_extra_array(AV);
}
