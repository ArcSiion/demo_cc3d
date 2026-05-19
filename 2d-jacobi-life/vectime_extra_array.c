#include "define.h"

void vectime_extra_array(int* A, int NX, int NY, int T) {

	int (* B)[NX + 2 * XSTART][ NY + 2 * YSTART] =  (int (*)[NX + 2 * XSTART][ NY + 2 * YSTART]) A;
    int tmp[8];

    int tt, t, x, xx, y;

    __m256i v_center_0, v_center_1, v_center_2; 
    __m256i v_x_plus_0, v_x_plus_1, v_x_plus_2; 
    __m256i v_x_minus_0, v_x_minus_1, v_x_minus_2;
    __m256i in, out, vzero, vone, vtwo, vthree, vrotatei_high2low;

    int zero[VECLEN_INT] =  {0, 0, 0, 0, 0, 0, 0, 0};
    int one[VECLEN_INT] =   {1, 1, 1, 1, 1, 1, 1, 1};
    int two[VECLEN_INT] =   {2, 2, 2, 2, 2, 2, 2, 2};
    int three[VECLEN_INT] = {3, 3, 3, 3, 3, 3, 3, 3};
 	int rotatei_high2low[VECLEN_INT] = {7, 0, 1, 2, 3, 4, 5, 6};



	int (* AV) [NY][VECLEN_INT] = (int(*)[NY][VECLEN_INT])alloc_extra_array(sizeof(int) * NY * VECLEN_INT * 3);

	int (* BV0) [VECLEN_INT] = (int (*) [VECLEN_INT]) AV;
	int (* BV1) [VECLEN_INT] = (int (*) [VECLEN_INT]) (AV + 1);
    int (* BV2) [VECLEN_INT] = (int (*) [VECLEN_INT]) (AV + 2);

    int (* Btmp [3]) [VECLEN_INT]  = {BV0, BV1, BV2};

	for ( tt = 0; tt <= T - VECLEN_INT; tt += VECLEN_INT){	

		for( t = tt ; t < tt + VECLEN_INT - 1 ; t++){
			for ( x = XSTART; x < XSTART + STRIDE * (VECLEN_INT - 1 - (t - tt)); x++) {
            #pragma ivdep
            #pragma vector always
                for ( y = YSTART; y < NY + YSTART; y++) {
                    Compute_scalar(B, t, x, y);
                }		
            }
		}

        t = tt;

        for(x = XSTART - XSLOPE; x <= XSTART +XSLOPE; x++){
			for ( y = YSTART; y <= NY + YSTART - VECLEN_INT; y += VECLEN_INT ){
				vloadi2(v_center_0, B[t%2    ][x + STRIDE * 7][y]);
				vloadi2(v_center_1, B[(t+1)%2][x + STRIDE * 6][y]);
				vloadi2(v_center_2, B[t%2    ][x + STRIDE * 5][y]);
				vloadi2(vrotatei_high2low, B[(t+1)%2][x + STRIDE * 4][y]);
				vloadi2(v_x_plus_0, B[t%2    ][x + STRIDE * 3][y]);
				vloadi2(v_x_plus_1, B[(t+1)%2][x + STRIDE * 2][y]);
				vloadi2(v_x_plus_2, B[t%2    ][x + STRIDE * 1][y]);
				vloadi2(vtwo, B[(t+1)%2][x             ][y]);
                transposei(v_center_0, v_center_1, v_center_2, vrotatei_high2low, v_x_plus_0, v_x_plus_1, v_x_plus_2, vtwo, in, out, vzero, vone, v_x_minus_0, v_x_minus_1, v_x_minus_2, vthree);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 7][0], vthree);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 6][0], v_x_minus_2);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 5][0], v_x_minus_1);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 4][0], v_x_minus_0);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 3][0], vone);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 2][0], vzero);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 1][0], out);
				vstorei( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 0][0], in);
			}
		}        

        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VECLEN_INT ; x ++){

            y = YSTART;
     	
	    	vzero = vloadi(zero[0]); 
            vone = vloadi(one[0]); 
            vtwo = vloadi(two[0]); 
            vthree = vloadi(three[0]);
            vrotatei_high2low = vloadi(rotatei_high2low[0]);

            vloadseti(v_x_minus_0, B, t, x-XSLOPE, y-YSLOPE);  
            vloadseti(v_center_0, B, t, x, y-YSLOPE);  
            vloadseti(v_x_plus_0, B, t, x+XSLOPE, y-YSLOPE); 


            vloadi2(v_x_minus_1, BV0[0][0]);
            vloadi2(v_center_1,  BV1[0][0]);
            vloadi2(v_x_plus_1,  BV2[0][0]);

            for ( y = YSTART; y <= NY + YSTART - VECLEN_INT; y+= VECLEN_INT) {

 
                vloadi2(in, B[(t)%2][x+STRIDE*VECLEN_INT][y]); 

                vloadi2(v_x_minus_2, BV0[y - YSTART + 1][0]);
                vloadi2(v_center_2,  BV1[y - YSTART + 1][0]); 
                vloadi2(v_x_plus_2,  BV2[y - YSTART + 1][0]); 
                Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_1(out, v_center_0, in);
                vstorei(BV0[y - YSTART][0], v_center_0); 


                vloadi2(v_x_minus_0, BV0[y - YSTART + 2][0]);
                vloadi2(v_center_0,  BV1[y - YSTART + 2][0]);
                vloadi2(v_x_plus_0,  BV2[y - YSTART + 2][0]);
                Compute_1vector(v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_2(out, v_center_1, in);
                vstorei(BV0[y - YSTART + 1][0], v_center_1); 


                vloadi2(v_x_minus_1, BV0[y - YSTART + 3][0]);
                vloadi2(v_center_1,  BV1[y - YSTART + 3][0]);
                vloadi2(v_x_plus_1,  BV2[y - YSTART + 3][0]);
                Compute_1vector(v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_3(out, v_center_2, in);	
                vstorei(BV0[y - YSTART + 2][0], v_center_2); 

                vloadi2(v_x_minus_2, BV0[y - YSTART + 4][0]);
                vloadi2(v_center_2,  BV1[y - YSTART + 4][0]);
                vloadi2(v_x_plus_2,  BV2[y - YSTART + 4][0]);
                Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_4(out, v_center_0, in);
                vstorei(BV0[y - YSTART + 3][0], v_center_0); 

                vloadi2(v_x_minus_0, BV0[y - YSTART + 5][0]);
                vloadi2(v_center_0,  BV1[y - YSTART + 5][0]);
                vloadi2(v_x_plus_0,  BV2[y - YSTART + 5][0]);
                Compute_1vector(v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_5(out, v_center_1, in);
                vstorei(BV0[y - YSTART + 4][0], v_center_1); 

                vloadi2(v_x_minus_1, BV0[y - YSTART + 6][0]);
                vloadi2(v_center_1,  BV1[y - YSTART + 6][0]);
                vloadi2(v_x_plus_1,  BV2[y - YSTART + 6][0]);
                Compute_1vector(v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_6(out, v_center_2, in);	
                vstorei(BV0[y - YSTART + 5][0], v_center_2); 


                vloadi2(v_x_minus_2, BV0[y - YSTART + 7][0]);
                vloadi2(v_center_2,  BV1[y - YSTART + 7][0]);
                vloadi2(v_x_plus_2,  BV2[y - YSTART + 7][0]);
                Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_7(out, v_center_0, in);
                vstorei(BV0[y - YSTART + 6][0], v_center_0); 



                if(y + VECLEN_INT <= NY + YSTART - VECLEN_INT){

                    vloadi2(v_x_minus_0, BV0[y - YSTART + 8][0]);
                    vloadi2(v_center_0,  BV1[y - YSTART + 8][0]);
                    vloadi2(v_x_plus_0,  BV2[y - YSTART + 8][0]);
                } else {
                    vloadseti(v_x_minus_0, B, t, x-XSLOPE, y+VECLEN_INT);  
                    vloadseti(v_center_0, B, t, x, y+VECLEN_INT);  
                    vloadseti(v_x_plus_0, B, t, x+XSLOPE, y+VECLEN_INT); 
                }
                Compute_1vector(v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_8(out, v_center_1, in);
                vstorei(BV0[y - YSTART + 7][0], v_center_1);

                v_x_minus_1 = v_x_minus_0;
                v_center_1 = v_center_0;
                v_x_plus_1 = v_x_plus_0;
                v_x_minus_0 = v_x_minus_2;
                v_center_0 = v_center_2;
                v_x_plus_0 = v_x_plus_2;


                vstorei(B[(t)%2][x][y], out);

            }
         
            for ( y += 1; y <= NY+YSTART; y++){
                vloadseti(v_x_minus_2, B, t, x-XSLOPE, y);  
                vloadseti(v_center_2, B, t, x, y);  
                vloadseti(v_x_plus_2, B, t, x+XSLOPE, y); 

                Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0); 
                vstorei(tmp[0], v_center_0);


		    	B[(t+1)%2][x + STRIDE * 7][y - 1] = tmp[0];
		    	B[t%2    ][x + STRIDE * 6][y - 1] = tmp[1];
		    	B[(t+1)%2][x + STRIDE * 5][y - 1] = tmp[2];
		    	B[t%2    ][x + STRIDE * 4][y - 1] = tmp[3];
		    	B[(t+1)%2][x + STRIDE * 3][y - 1] = tmp[4];
		    	B[t%2    ][x + STRIDE * 2][y - 1] = tmp[5];
		    	B[(t+1)%2][x + STRIDE * 1][y - 1] = tmp[6];
		    	B[t%2    ][x             ][y - 1] = tmp[7];

                v_x_minus_0 = v_x_minus_1;
                v_center_0 = v_center_1;
                v_x_plus_0 = v_x_plus_1;
                v_x_minus_1 = v_x_minus_2;
                v_center_1 = v_center_2;
                v_x_plus_1 = v_x_plus_2;
                
            }
            Btmp[0] = BV0;
            BV0 = BV1;
            BV1 = BV2;
            BV2 = Btmp[0]; 
        }
        Btmp [0] = BV0;
        Btmp [1] = BV1;
        Btmp [2] = BV2;

        for ( y = YSTART; y < NY + YSTART ; y++ ){
            B[1][XSTART - XSLOPE][y] = B[0][XSTART - XSLOPE][y];
        }

        xx = x;
        for(x = xx - XSLOPE; x <= xx +XSLOPE; x++){
			for ( y = YSTART; y <= NY + YSTART - VECLEN_INT; y += VECLEN_INT ){

				vloadi2(vtwo, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 7][0]);
				vloadi2(v_x_plus_2, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 6][0]);
				vloadi2(v_x_plus_1, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 5][0]);
				vloadi2(v_x_plus_0, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 4][0]);
				vloadi2(vrotatei_high2low, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 3][0]);
				vloadi2(v_center_2, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 2][0]);
				vloadi2(v_center_1, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 1][0]);
				vloadi2(v_center_0, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 0][0]);
                transposei(v_center_0, v_center_1, v_center_2, vrotatei_high2low, v_x_plus_0, v_x_plus_1, v_x_plus_2, vtwo, in, out, vzero, vone, v_x_minus_0, v_x_minus_1, v_x_minus_2, vthree);
				vstorei( B[t%2    ][x + STRIDE *7][y], in);
				vstorei( B[(t+1)%2][x + STRIDE *6][y], out);
				vstorei( B[t%2    ][x + STRIDE *5][y], vzero);
				vstorei( B[(t+1)%2][x + STRIDE *4][y], vone);
				vstorei( B[t%2    ][x + STRIDE *3][y], v_x_minus_0);
				vstorei( B[(t+1)%2][x + STRIDE *2][y], v_x_minus_1);
				vstorei( B[t%2    ][x + STRIDE *1][y], v_x_minus_2);
				vstorei( B[(t+1)%2][x            ][y], vthree);
			}
		} 
    

        for( t = tt ; t < tt + VECLEN_INT ; t++){
			for ( x = xx + STRIDE * (VECLEN_INT - 1 - (t - tt)); x < NX + XSTART; x++) {
                #pragma ivdep
                #pragma vector always
                for ( y = YSTART; y < NY + YSTART; y++) {
                    Compute_scalar(B, t, x, y);
                }		
            }
		}
	}

	for (t = tt ; t < T; t++){
		for (x = XSTART; x < NX + XSTART; x++) {
            #pragma ivdep
            #pragma vector always
            for ( y = YSTART; y < NY + YSTART; y++) {
               Compute_scalar(B, t, x, y);
            }
		}
	}	
    free_extra_array(AV);
}
