#include "define_512.h"

void vectime_extra_array_512(double* A, int NX, int NY, int NZ, int T) {
	if (back_scalar_512(NX, NY, NZ, T)) {
		naive_scalar(A, NX, NY, NZ, T);
		return;
	}
	
    double (* B)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;

    long int i, j, k;

    double tmp[VVECLEN];
    int tt, t = 0, x, xx, y, yy, z, zz;

    vec_512 v_y_plus_1, v_y_minus_1;
    vec_512 v_center_0, v_center_1, v_center_2, v_center_3, v_center_4;
    vec_512 v_center_5, v_center_6, v_center_7, v_center_8, v_center_9;

    vec_512 in, out;
    SET_COFF_512;

    double (* AV) [NY][NZ][VVECLEN] = (double(*)[NY][NZ][VVECLEN])alloc_extra_array_512(sizeof(double) * NY * NZ * VVECLEN * 3);

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
                    for ( z = ZSTART; z < NZ+ZSTART; z++) {
                        Compute_scalar(B, t, x, y, z);
                    }   
                }       
            }
        }
        t = tt;

        // 转置
        for(x = XSTART - XSLOPE; x <= XSTART + XSLOPE; x++){
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
                    
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 0][0], v_center_0);
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 1][0], v_center_1);
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 2][0], v_center_2);
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 3][0], v_center_3);
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 4][0], v_center_4);
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 5][0], v_center_5);
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 6][0], v_center_6);
                    vstore_512(Btmp[x - XSTART + XSLOPE][y - YSTART][(z - ZSTART) + 7][0], v_center_7);
                }
            }
        } 
        
        // 主体计算推进
        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VVECLEN ; x ++){

            for ( y = YSTART; y < NY + YSTART; y++ ) {
                z = ZSTART;
                
                vset_3d_512(v_center_0, x, y, z-ZSLOPE);  // z-1
                vload_512(v_center_1, BV1[y - YSTART][0][0]); // z

                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                
                    vload_512(in, B[(t)%2][x+STRIDE*VVECLEN][y][z]);   // the next x iter in vector

                    // =================================== 1st ===================================
                    vload_512(v_center_2 , BV1[y - YSTART][z - ZSTART + 1][0]);                                                                        
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART][0]));   // (0,-1,0)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z)   \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART][0]));   // (0,+1,0)
                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART][0])));   // 1st   
                    Input_Output_1_512(out, v_center_0, in);
                    vstore_512(BV0[y - YSTART][z - ZSTART][0], v_center_0);
                    
                    // =================================== 2nd ===================================
                    vload_512(v_center_3 , BV1[y - YSTART][z - ZSTART + 2][0]);                  // (0,0,+2)                                                   
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z + 1)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART + 1][0]));    // (0,-1,+1)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z + 1)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART + 1][0]));    // (0,+1,+1)          
                    Compute_1vector_512(v_center_2, v_center_1, v_center_3, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART + 1][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART + 1][0])));   // 2nd   
                    Input_Output_2_512(out, v_center_1, in);
                    vstore_512(BV0[y - YSTART][z - ZSTART + 1][0], v_center_1); 
                    
                    // =================================== 3rd ===================================
                    vload_512(v_center_4 , BV1[y - YSTART][z - ZSTART + 3][0]);                    // (0,0,+3)                                                   
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z + 2)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART + 2][0]));    // (0,-1,+2)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z + 2)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART + 2][0]));    // (0,+1,+2)
                    Compute_1vector_512(v_center_3, v_center_2, v_center_4, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART + 2][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART + 2][0])));   // 3rd   
                    Input_Output_3_512(out, v_center_2, in);    
                    vstore_512(BV0[y - YSTART][z - ZSTART + 2][0], v_center_2); 
                    
                    // =================================== 4th ===================================
                    vload_512(v_center_5 , BV1[y - YSTART][z - ZSTART + 4][0]);                    // (0,0,+4)                                                   
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z + 3)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART + 3][0]));    // (0,-1,+3)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z + 3)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART + 3][0]));    // (0,+1,+3)
                    Compute_1vector_512(v_center_4, v_center_3, v_center_5, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART + 3][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART + 3][0])));   // 4th   
                    Input_Output_4_512(out, v_center_3, in);    
                    vstore_512(BV0[y - YSTART][z - ZSTART + 3][0], v_center_3); 
                    
                    // =================================== 5th ===================================
                    vload_512(v_center_6 , BV1[y - YSTART][z - ZSTART + 5][0]);                    // (0,0,+5)                                                   
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z + 4)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART + 4][0]));    
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z + 4)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART + 4][0]));    
                    Compute_1vector_512(v_center_5, v_center_4, v_center_6, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART + 4][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART + 4][0])));   // 5th   
                    Input_Output_5_512(out, v_center_4, in);    
                    vstore_512(BV0[y - YSTART][z - ZSTART + 4][0], v_center_4); 

                    // =================================== 6th ===================================
                    vload_512(v_center_7 , BV1[y - YSTART][z - ZSTART + 6][0]);                    // (0,0,+6)                                                   
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z + 5)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART + 5][0]));    
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z + 5)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART + 5][0]));    
                    Compute_1vector_512(v_center_6, v_center_5, v_center_7, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART + 5][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART + 5][0])));   // 6th   
                    Input_Output_6_512(out, v_center_5, in);    
                    vstore_512(BV0[y - YSTART][z - ZSTART + 5][0], v_center_5); 

                    // =================================== 7th ===================================
                    vload_512(v_center_8 , BV1[y - YSTART][z - ZSTART + 7][0]);                    // (0,0,+7)                                                   
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z + 6)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART + 6][0]));    
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z + 6)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART + 6][0]));    
                    Compute_1vector_512(v_center_7, v_center_6, v_center_8, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART + 6][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART + 6][0])));   // 7th   
                    Input_Output_7_512(out, v_center_6, in);    
                    vstore_512(BV0[y - YSTART][z - ZSTART + 6][0], v_center_6); 

                    // =================================== 8th ===================================
                    v_center_9 = ( z > NZ + ZSTART - VVECLEN - VVECLEN ) ? \
                                    setv_3d_512(x, y, z + 8) \
                                :   _mm512_loadu_pd(&(BV1[y - YSTART][z - ZSTART + 8][0]))  ;               // (0,0,+8) 
                    v_y_minus_1 = ( y == YSTART ) ? \
                                        setv_3d_512(x, y - YSLOPE, z + 7)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART - YSLOPE][z - ZSTART + 7][0]));    // (0,-1,+7)
                    v_y_plus_1  = ( y == NY + YSTART - YSLOPE ) ? \
                                        setv_3d_512(x, y + YSLOPE, z + 7)    \
                                    :   _mm512_loadu_pd(&(BV1[y - YSTART + YSLOPE][z - ZSTART + 7][0]));    // (0,+1,+7)
                    Compute_1vector_512(v_center_8, v_center_7, v_center_9, v_y_minus_1, \
                                    v_y_plus_1, \
                                    _mm512_loadu_pd(&(BV0[y - YSTART][z - ZSTART + 7][0])), \
                                    _mm512_loadu_pd(&(BV2[y - YSTART][z - ZSTART + 7][0])));   // 8th   

                    Input_Output_8_512(out, v_center_7, in);    
                    vstore_512(BV0[y - YSTART][z - ZSTART + 7][0], v_center_7); 
                    //-------------------------------------------------------------------------------------------

                    // out 包含了 x 切片上的结果，回写到主存中
                    vstore_512(B[(t)%2][x][y][z], out);
                    
                    v_center_0 = v_center_8;
                    v_center_1 = v_center_9;
                }
                
                // 处理残局标量 Z
                for ( ; z < NZ + ZSTART; z++) {
                    vset_3d_512(v_center_2 , x, y, z + 1);

                    Compute_1vector_512(v_center_1, v_center_0, v_center_2, \
                                    setv_3d_512(x, y - 1, z), \
                                    setv_3d_512(x, y + 1, z), \
                                    setv_3d_512(x - 1, y, z), \
                                    setv_3d_512(x + 1, y, z));
                    
                    vstore_set_3d_512(x, y, z, v_center_0);

                    v_center_0 = v_center_1;
                    v_center_1 = v_center_2;
                } 
            }
           
            // 只有在 X 层推进一步时，才对 BV 数组进行切片替换轮转
            Btmp[0] = BV0;
            BV0 = BV1;
            BV1 = BV2;
            BV2 = Btmp[0];     
        }
        
        // 恢复数组状态为这轮计算的结尾形态
        Btmp [0] = BV0;
        Btmp [1] = BV1;
        Btmp [2] = BV2;
        
        // 将缓存在 Btmp 里面的结果写回正确的倾斜内存区
        for(; x < NX + XSTART - STRIDE * VVECLEN + 1 + 3; x++){           
            for ( y = YSTART ; y < NY + YSTART ; y ++ ){
                for ( z = ZSTART; z <= NZ + ZSTART - VVECLEN; z += VVECLEN) {
                    vload_512(v_center_0, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 0][0]);
                    vload_512(v_center_1, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 1][0]);
                    vload_512(v_center_2, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 2][0]);
                    vload_512(v_center_3, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 3][0]);    
                    vload_512(v_center_4, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 4][0]);
                    vload_512(v_center_5, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 5][0]);
                    vload_512(v_center_6, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 6][0]);
                    vload_512(v_center_7, Btmp[x - (NX + XSTART - STRIDE * VVECLEN + 1)][(y - YSTART)][(z - ZSTART) + 7][0]);
                    
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

        // Tail
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
    
    // Extra points 处理末尾时间步
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
    
    // 释放额外数组内存
    free_extra_array(AV);
}
