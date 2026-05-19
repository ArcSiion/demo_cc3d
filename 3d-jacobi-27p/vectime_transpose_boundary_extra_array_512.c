#include "define_512.h"

void vectime_transpose_boundary_extra_array_512(double* A, int NX, int NY, int NZ, int T) {
    if (back_scalar_512(NX, NY, NZ, T)) {
        naive_scalar(A, NX, NY, NZ, T);
        return;
    }
    
    double (* B)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;

    long int i, j, k;

    double tmp[VVECLEN];
    int tt, t = 0, x, xx, y, yy, z, zz;
    
    vec_512 v_center_0, v_center_1, v_center_2, v_center_3, v_center_4;
    vec_512 v_center_5, v_center_6, v_center_7, v_center_8, v_center_9;
    
    vec_512 v_all_d_1_0, v_all_d_1_1, v_all_d_1_2, v_all_d_1_3, v_all_d_1_4;
    vec_512 v_all_d_1_5, v_all_d_1_6, v_all_d_1_7, v_all_d_1_8, v_all_d_1_9;
    
    vec_512 v_all_d_2_0, v_all_d_2_1, v_all_d_2_2, v_all_d_2_3, v_all_d_2_4;
    vec_512 v_all_d_2_5, v_all_d_2_6, v_all_d_2_7, v_all_d_2_8, v_all_d_2_9;

    vec_512 in, out;
    SET_COFF_512;

    double (* AV) [NY + 2][NZ][VVECLEN] = (double(*)[NY+2][NZ][VVECLEN])alloc_extra_array_512(sizeof(double) * (NY + 2) * NZ * VVECLEN * 4);

    double (* BV0) [NZ][VVECLEN] = (double (*) [NZ][VVECLEN]) AV;
    double (* BV1) [NZ][VVECLEN] = (double (*) [NZ][VVECLEN]) (AV + 1);
    double (* BV2) [NZ][VVECLEN] = (double (*) [NZ][VVECLEN]) (AV + 2);
    double (* BV3) [NZ][VVECLEN] = (double (*) [NZ][VVECLEN]) (AV + 3);

    double (* Btmp [4]) [NZ][VVECLEN]  = {BV0, BV1, BV2, BV3};

    for ( tt = 0; tt <= T - VVECLEN; tt += VVECLEN){  
        for( t = tt ; t < tt + VVECLEN - 1 ; t++){      //head
            for ( x = XSTART; x < XSTART + STRIDE * (VVECLEN - 1 - (t - tt)); x++) { //ASSERT VVECLEN <= STRIDE + 1
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

        // 初始数据转置 8x8
        for(x = XSTART - XSLOPE; x <= XSTART + XSLOPE; x++){
            for ( y = YSTART - YSLOPE; y <= NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ; z += VVECLEN ) {
                    vload_512(v_center_0, B[t%2    ][x + STRIDE *7][y][z]);
                    vload_512(v_center_1, B[(t+1)%2][x + STRIDE *6][y][z]);
                    vload_512(v_center_2, B[t%2    ][x + STRIDE *5][y][z]);
                    vload_512(v_center_3, B[(t+1)%2][x + STRIDE *4][y][z]);
                    vload_512(v_center_4, B[t%2    ][x + STRIDE *3][y][z]);
                    vload_512(v_center_5, B[(t+1)%2][x + STRIDE *2][y][z]);
                    vload_512(v_center_6, B[t%2    ][x + STRIDE *1][y][z]);
                    vload_512(v_center_7, B[(t+1)%2][x            ][y][z]);
                    
                    transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);
                    
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 0][0], v_center_0);
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 1][0], v_center_1);
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 2][0], v_center_2);
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 3][0], v_center_3);
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 4][0], v_center_4);
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 5][0], v_center_5);
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 6][0], v_center_6);
                    _mm512_storeu_pd(&Btmp[x - XSTART + XSLOPE][y - YSTART + YSLOPE][(z - ZSTART) + 7][0], v_center_7);
                }
            }
        }

        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VVECLEN ; x ++){
 
            for ( y = YSTART; y < NY + YSTART; y++ ) {

                z = ZSTART;

                // Loop 启动前的 Z-1 层初始化
                v_center_0  = load_v2_512(x, y, z, 0, 0, -1);
                v_all_d_1_0 = Add_4_vectors_512(load_v2_512(x, y, z, 1, 0, -1), load_v2_512(x, y, z, 0, 1, -1), load_v2_512(x, y, z, -1, 0, -1), load_v2_512(x, y, z, 0, -1, -1));
                v_all_d_2_0 = Add_4_vectors_512(load_v2_512(x, y, z, 1, 1, -1), load_v2_512(x, y, z, -1, 1, -1), load_v2_512(x, y, z, -1, -1, -1), load_v2_512(x, y, z, 1, -1, -1));

                // Loop 启动前的 Z 层初始化
                v_center_1  = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][0][0]);
                v_all_d_1_1 = Add_4_vectors_512(load_x_m_512(x, y, z, 0, 0), load_x_p_512(x, y, z, 0, 0), load_x_c_512(x, y, z, 1, 0), load_x_c_512(x, y, z, -1, 0));
                v_all_d_2_1 = Add_4_vectors_512(load_x_m_512(x, y, z, 1, 0), load_x_m_512(x, y, z, -1, 0), load_x_p_512(x, y, z, 1, 0), load_x_p_512(x, y, z, -1, 0));

                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                
                    vload_512(in, B[(t)%2][x+STRIDE*VVECLEN][y][z]);   

                    // ------------------ Step 1 (z+0) ------------------
                    v_center_2 = load_x_c_512(x, y, z, 0, 1); 
                    v_all_d_1_2 = Add_4_vectors_512(load_x_m_512(x,y,z,0,1), load_x_p_512(x,y,z,0,1), load_x_c_512(x,y,z,1,1), load_x_c_512(x,y,z,-1,1));
                    v_all_d_2_2 = Add_4_vectors_512(load_x_m_512(x,y,z,1,1), load_x_m_512(x,y,z,-1,1), load_x_p_512(x,y,z,1,1), load_x_p_512(x,y,z,-1,1));
                    
                    Compute_1vector_512(v_center_0, v_center_1, v_center_2, v_all_d_1_0, v_all_d_1_1, v_all_d_1_2, v_all_d_2_0, v_all_d_2_1, v_all_d_2_2);
                    Input_Output_1_512(out, v_center_0, in);
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART][0], v_center_0);
                    
                    // ------------------ Step 2 (z+1) ------------------
                    v_center_3 = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 2][0]);
                    v_all_d_1_3 = Add_4_vectors_512(load_x_m_512(x,y,z,0,2), load_x_p_512(x,y,z,0,2), load_x_c_512(x,y,z,1,2), load_x_c_512(x,y,z,-1,2));
                    v_all_d_2_3 = Add_4_vectors_512(load_x_m_512(x,y,z,1,2), load_x_m_512(x,y,z,-1,2), load_x_p_512(x,y,z,1,2), load_x_p_512(x,y,z,-1,2));
                    
                    Compute_1vector_512(v_center_1, v_center_2, v_center_3, v_all_d_1_1, v_all_d_1_2, v_all_d_1_3, v_all_d_2_1, v_all_d_2_2, v_all_d_2_3);
                    Input_Output_2_512(out, v_center_1, in);
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART + 1][0], v_center_1);
                    
                    // ------------------ Step 3 (z+2) ------------------
                    v_center_4 = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 3][0]);
                    v_all_d_1_4 = Add_4_vectors_512(load_x_m_512(x,y,z,0,3), load_x_p_512(x,y,z,0,3), load_x_c_512(x,y,z,1,3), load_x_c_512(x,y,z,-1,3));
                    v_all_d_2_4 = Add_4_vectors_512(load_x_m_512(x,y,z,1,3), load_x_m_512(x,y,z,-1,3), load_x_p_512(x,y,z,1,3), load_x_p_512(x,y,z,-1,3));
                    
                    Compute_1vector_512(v_center_2, v_center_3, v_center_4, v_all_d_1_2, v_all_d_1_3, v_all_d_1_4, v_all_d_2_2, v_all_d_2_3, v_all_d_2_4);
                    Input_Output_3_512(out, v_center_2, in);
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART + 2][0], v_center_2);
                    
                    // ------------------ Step 4 (z+3) ------------------
                    v_center_5 = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 4][0]);
                    v_all_d_1_5 = Add_4_vectors_512(load_x_m_512(x,y,z,0,4), load_x_p_512(x,y,z,0,4), load_x_c_512(x,y,z,1,4), load_x_c_512(x,y,z,-1,4));
                    v_all_d_2_5 = Add_4_vectors_512(load_x_m_512(x,y,z,1,4), load_x_m_512(x,y,z,-1,4), load_x_p_512(x,y,z,1,4), load_x_p_512(x,y,z,-1,4));
                    
                    Compute_1vector_512(v_center_3, v_center_4, v_center_5, v_all_d_1_3, v_all_d_1_4, v_all_d_1_5, v_all_d_2_3, v_all_d_2_4, v_all_d_2_5);
                    Input_Output_4_512(out, v_center_3, in);
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART + 3][0], v_center_3);
                    
                    // ------------------ Step 5 (z+4) ------------------
                    v_center_6 = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 5][0]);
                    v_all_d_1_6 = Add_4_vectors_512(load_x_m_512(x,y,z,0,5), load_x_p_512(x,y,z,0,5), load_x_c_512(x,y,z,1,5), load_x_c_512(x,y,z,-1,5));
                    v_all_d_2_6 = Add_4_vectors_512(load_x_m_512(x,y,z,1,5), load_x_m_512(x,y,z,-1,5), load_x_p_512(x,y,z,1,5), load_x_p_512(x,y,z,-1,5));
                    
                    Compute_1vector_512(v_center_4, v_center_5, v_center_6, v_all_d_1_4, v_all_d_1_5, v_all_d_1_6, v_all_d_2_4, v_all_d_2_5, v_all_d_2_6);
                    Input_Output_5_512(out, v_center_4, in);
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART + 4][0], v_center_4);

                    // ------------------ Step 6 (z+5) ------------------
                    v_center_7 = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 6][0]);
                    v_all_d_1_7 = Add_4_vectors_512(load_x_m_512(x,y,z,0,6), load_x_p_512(x,y,z,0,6), load_x_c_512(x,y,z,1,6), load_x_c_512(x,y,z,-1,6));
                    v_all_d_2_7 = Add_4_vectors_512(load_x_m_512(x,y,z,1,6), load_x_m_512(x,y,z,-1,6), load_x_p_512(x,y,z,1,6), load_x_p_512(x,y,z,-1,6));
                    
                    Compute_1vector_512(v_center_5, v_center_6, v_center_7, v_all_d_1_5, v_all_d_1_6, v_all_d_1_7, v_all_d_2_5, v_all_d_2_6, v_all_d_2_7);
                    Input_Output_6_512(out, v_center_5, in);
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART + 5][0], v_center_5);

                    // ------------------ Step 7 (z+6) ------------------
                    v_center_8 = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 7][0]);
                    v_all_d_1_8 = Add_4_vectors_512(load_x_m_512(x,y,z,0,7), load_x_p_512(x,y,z,0,7), load_x_c_512(x,y,z,1,7), load_x_c_512(x,y,z,-1,7));
                    v_all_d_2_8 = Add_4_vectors_512(load_x_m_512(x,y,z,1,7), load_x_m_512(x,y,z,-1,7), load_x_p_512(x,y,z,1,7), load_x_p_512(x,y,z,-1,7));
                    
                    Compute_1vector_512(v_center_6, v_center_7, v_center_8, v_all_d_1_6, v_all_d_1_7, v_all_d_1_8, v_all_d_2_6, v_all_d_2_7, v_all_d_2_8);
                    Input_Output_7_512(out, v_center_6, in);
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART + 6][0], v_center_6);

                    // ------------------ Step 8 (z+7) ------------------
                    // 处理跨边界加载问题
                    if ( z > NZ + ZSTART - VVECLEN - VVECLEN ){
                        v_center_9  = load_v2_512(x, y, z, 0, 0, VVECLEN);
                        v_all_d_1_9 = Add_4_vectors_512(load_v2_512(x,y,z,1,0,VVECLEN), load_v2_512(x,y,z,0,1,VVECLEN), load_v2_512(x,y,z,-1,0,VVECLEN), load_v2_512(x,y,z,0,-1,VVECLEN));
                        v_all_d_2_9 = Add_4_vectors_512(load_v2_512(x,y,z,1,1,VVECLEN), load_v2_512(x,y,z,-1,1,VVECLEN), load_v2_512(x,y,z,-1,-1,VVECLEN), load_v2_512(x,y,z,1,-1,VVECLEN));
                    } else {
                        v_center_9  = _mm512_loadu_pd( & BV1[y - YSTART + YSLOPE][z - ZSTART + 8][0]);
                        v_all_d_1_9 = Add_4_vectors_512(load_x_m_512(x,y,z,0,8), load_x_p_512(x,y,z,0,8), load_x_c_512(x,y,z,1,8), load_x_c_512(x,y,z,-1,8));
                        v_all_d_2_9 = Add_4_vectors_512(load_x_m_512(x,y,z,1,8), load_x_m_512(x,y,z,-1,8), load_x_p_512(x,y,z,1,8), load_x_p_512(x,y,z,-1,8));
                    }

                    Compute_1vector_512(v_center_7, v_center_8, v_center_9, v_all_d_1_7, v_all_d_1_8, v_all_d_1_9, v_all_d_2_7, v_all_d_2_8, v_all_d_2_9);
                    Input_Output_8_512(out, v_center_7, in);    
                    vstore_512(BV3[y - YSTART + YSLOPE][z - ZSTART + 7][0], v_center_7);
                    
                    vstore_512(B[(t)%2][x][y][z], out);

                    // 寄存器移位，为下一个 Z-Block 初始化
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                    
                    v_all_d_1_0 = v_all_d_1_8;
                    v_all_d_1_1 = v_all_d_1_9;
                    
                    v_all_d_2_0 = v_all_d_2_8;
                    v_all_d_2_1 = v_all_d_2_9;
                }
                for ( ; z < NZ + ZSTART; z++) {

                    v_center_2 = setv_3d_512(x, y, z + 1);
                    v_all_d_1_2 = Add_4_vectors_512(
                        setv_3d_512(x + 1, y    , z + 1),
                        setv_3d_512(x    , y + 1, z + 1),
                        setv_3d_512(x - 1, y    , z + 1),
                        setv_3d_512(x    , y - 1, z + 1)
                    );
                    v_all_d_2_2 = Add_4_vectors_512(
                        setv_3d_512(x + 1, y + 1, z + 1),
                        setv_3d_512(x - 1, y + 1, z + 1),
                        setv_3d_512(x - 1, y - 1, z + 1),
                        setv_3d_512(x + 1, y - 1, z + 1)
                    );

                    Compute_1vector_512(
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

                    vstore_set_3d_512(x, y, z, v_center_0);

                    v_center_0  = v_center_1;
                    v_center_1  = v_center_2;
                    v_all_d_1_0 = v_all_d_1_1;
                    v_all_d_1_1 = v_all_d_1_2;
                    v_all_d_2_0 = v_all_d_2_1;
                    v_all_d_2_1 = v_all_d_2_2;
                }
            }

            // 处理 Y 轴下边界缓冲的转置写入
            y = YSTART - YSLOPE;
            for ( z = ZSTART ; z <= NZ + ZSTART - VVECLEN; z += VVECLEN){ 
                vload_512(v_center_0, B[t%2    ][x + XSLOPE + 1 + STRIDE *7][YSTART - YSLOPE][z]);
                vload_512(v_center_1, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *6][YSTART - YSLOPE][z]);
                vload_512(v_center_2, B[t%2    ][x + XSLOPE + 1 + STRIDE *5][YSTART - YSLOPE][z]);
                vload_512(v_center_3, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *4][YSTART - YSLOPE][z]);
                vload_512(v_center_4, B[t%2    ][x + XSLOPE + 1 + STRIDE *3][YSTART - YSLOPE][z]);
                vload_512(v_center_5, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][YSTART - YSLOPE][z]);
                vload_512(v_center_6, B[t%2    ][x + XSLOPE + 1 + STRIDE *1][YSTART - YSLOPE][z]);
                vload_512(v_center_7, B[(t+1)%2][x + XSLOPE + 1            ][YSTART - YSLOPE][z]);

                // 处理 Y 轴上边界缓冲
                vload_512(v_all_d_1_0, B[t%2    ][x + XSLOPE + 1 + STRIDE *7][YSTART + NY][z]);
                vload_512(v_all_d_1_1, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *6][YSTART + NY][z]);
                vload_512(v_all_d_1_2, B[t%2    ][x + XSLOPE + 1 + STRIDE *5][YSTART + NY][z]);
                vload_512(v_all_d_1_3, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *4][YSTART + NY][z]);
                vload_512(v_all_d_1_4, B[t%2    ][x + XSLOPE + 1 + STRIDE *3][YSTART + NY][z]);
                vload_512(v_all_d_1_5, B[(t+1)%2][x + XSLOPE + 1 + STRIDE *2][YSTART + NY][z]);
                vload_512(v_all_d_1_6, B[t%2    ][x + XSLOPE + 1 + STRIDE *1][YSTART + NY][z]);
                vload_512(v_all_d_1_7, B[(t+1)%2][x + XSLOPE + 1            ][YSTART + NY][z]);

                transpose_512(v_center_0, v_center_1, v_center_2, v_center_3, v_center_4, v_center_5, v_center_6, v_center_7, in, out);
                transpose_512(v_all_d_1_0, v_all_d_1_1, v_all_d_1_2, v_all_d_1_3, v_all_d_1_4, v_all_d_1_5, v_all_d_1_6, v_all_d_1_7, in, out);

                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 7][0], v_center_7);
                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 6][0], v_center_6);
                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 5][0], v_center_5);
                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 4][0], v_center_4);
                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 3][0], v_center_3);
                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 2][0], v_center_2);
                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 1][0], v_center_1);
                _mm512_storeu_pd(&BV3[0][(z - ZSTART) + 0][0], v_center_0);

                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 0][0], v_all_d_1_0);
                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 1][0], v_all_d_1_1);
                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 2][0], v_all_d_1_2);
                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 3][0], v_all_d_1_3);
                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 4][0], v_all_d_1_4);
                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 5][0], v_all_d_1_5);
                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 6][0], v_all_d_1_6);
                _mm512_storeu_pd(&BV3[NY + YSLOPE][(z - ZSTART) + 7][0], v_all_d_1_7);
            }
            
            
            // 滚动数组切片轮转
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

        // 收尾区的 8x8 转置回写主存
        for(; x < NX + XSTART - STRIDE * VVECLEN + 1 + 3; x++){           
            for ( y = YSTART - YSLOPE ; y <= NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                    v_center_0 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 0][0]);
                    v_center_1 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 1][0]);
                    v_center_2 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 2][0]);
                    v_center_3 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 3][0]);  
                    v_center_4 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 4][0]);
                    v_center_5 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 5][0]);
                    v_center_6 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 6][0]);
                    v_center_7 = _mm512_loadu_pd(&Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART + YSLOPE)][(z - ZSTART) + 7][0]);  
                    
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
        
        //// tail 标量扫尾
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
    free_extra_array_512(AV);
}
