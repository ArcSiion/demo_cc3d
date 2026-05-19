#include "define_512.h"

void vectime_extra_array_512(int* A, int NX, int NY, int T) {

	int (* B)[NX + 2 * XSTART][ NY + 2 * YSTART] =  (int (*)[NX + 2 * XSTART][ NY + 2 * YSTART]) A;
    int tmp[VECLEN_INT_512];

    int tt, t, x, xx, y;

    veci_512 v_center_0, v_center_1, v_center_2; 
    veci_512 v_x_plus_0, v_x_plus_1, v_x_plus_2; 
    veci_512 v_x_minus_0, v_x_minus_1, v_x_minus_2;
    veci_512 in, out, vzero, vone, vtwo, vthree;

    veci_512 v0, v1, v2, v3, v4, v5, v6, v7;
    veci_512 v8, v9, v10, v11, v12, v13, v14, v15;
    SET_LIFE_CONST_512;


	int (* AV) [NY][VECLEN_INT_512] = (int(*)[NY][VECLEN_INT_512])alloc_extra_array_512(sizeof(int) * NY * VECLEN_INT_512 * 3);

	int (* BV0) [VECLEN_INT_512] = (int (*) [VECLEN_INT_512]) AV;
	int (* BV1) [VECLEN_INT_512] = (int (*) [VECLEN_INT_512]) (AV + 1);
    int (* BV2) [VECLEN_INT_512] = (int (*) [VECLEN_INT_512]) (AV + 2);

    int (* Btmp [3]) [VECLEN_INT_512]  = {BV0, BV1, BV2};

	for ( tt = 0; tt <= T - VECLEN_INT_512; tt += VECLEN_INT_512){	

		for( t = tt ; t < tt + VECLEN_INT_512 - 1 ; t++){
			for ( x = XSTART; x < XSTART + STRIDE * (VECLEN_INT_512 - 1 - (t - tt)); x++) {
            #pragma ivdep
            #pragma vector always
                for ( y = YSTART; y < NY + YSTART; y++) {
                    Compute_scalar(B, t, x, y);
                }		
            }
		}

        t = tt;

        for(x = XSTART - XSLOPE; x <= XSTART + XSLOPE; x++){
			for ( y = YSTART; y <= NY + YSTART - VECLEN_INT_512; y += VECLEN_INT_512 ){

				vloadi2_512(v0,  B[t%2    ][x + STRIDE * 15][y]);
				vloadi2_512(v1,  B[(t+1)%2][x + STRIDE * 14][y]);
				vloadi2_512(v2,  B[t%2    ][x + STRIDE * 13][y]);
				vloadi2_512(v3,  B[(t+1)%2][x + STRIDE * 12][y]);
				vloadi2_512(v4,  B[t%2    ][x + STRIDE * 11][y]);
				vloadi2_512(v5,  B[(t+1)%2][x + STRIDE * 10][y]);
				vloadi2_512(v6,  B[t%2    ][x + STRIDE * 9 ][y]);
				vloadi2_512(v7,  B[(t+1)%2][x + STRIDE * 8 ][y]);
				vloadi2_512(v8,  B[t%2    ][x + STRIDE * 7 ][y]);
				vloadi2_512(v9,  B[(t+1)%2][x + STRIDE * 6 ][y]);
				vloadi2_512(v10, B[t%2    ][x + STRIDE * 5 ][y]);
				vloadi2_512(v11, B[(t+1)%2][x + STRIDE * 4 ][y]);
				vloadi2_512(v12, B[t%2    ][x + STRIDE * 3 ][y]);
				vloadi2_512(v13, B[(t+1)%2][x + STRIDE * 2 ][y]);
				vloadi2_512(v14, B[t%2    ][x + STRIDE * 1 ][y]);
				vloadi2_512(v15, B[(t+1)%2][x             ][y]);

                transposei_512(v0, v1, v2, v3, v4, v5, v6, v7,
                               v8, v9, v10, v11, v12, v13, v14, v15,
                               v0, v1, v2, v3, v4, v5, v6, v7,
                               v8, v9, v10, v11, v12, v13, v14, v15);

				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 15][0], v15);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 14][0], v14);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 13][0], v13);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 12][0], v12);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 11][0], v11);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 10][0], v10);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 9 ][0], v9);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 8 ][0], v8);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 7 ][0], v7);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 6 ][0], v6);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 5 ][0], v5);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 4 ][0], v4);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 3 ][0], v3);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 2 ][0], v2);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 1 ][0], v1);
				vstorei_512( Btmp[x - (XSTART - XSLOPE)][(y - YSTART) + 0 ][0], v0);
			}
		}        

        for ( x = XSTART ; x <= NX + XSTART - STRIDE * VECLEN_INT_512 ; x ++){

            y = YSTART;
     	
            

            vloadseti_512(v_x_minus_0, B, t, x-XSLOPE, y-YSLOPE);  
            vloadseti_512(v_center_0, B, t, x, y-YSLOPE);  
            vloadseti_512(v_x_plus_0, B, t, x+XSLOPE, y-YSLOPE); 


            vloadi2_512(v_x_minus_1, BV0[0][0]);
            vloadi2_512(v_center_1,  BV1[0][0]);
            vloadi2_512(v_x_plus_1,  BV2[0][0]);

            for ( y = YSTART; y <= NY + YSTART - VECLEN_INT_512; y+= VECLEN_INT_512) {

 
                vloadi2_512(in, B[(t)%2][x+STRIDE*VECLEN_INT_512][y]); 

                vloadi2_512(v_x_minus_2, BV0[y - YSTART + 1][0]);
                vloadi2_512(v_center_2,  BV1[y - YSTART + 1][0]); 
                vloadi2_512(v_x_plus_2,  BV2[y - YSTART + 1][0]); 
                Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_1_512(out, v_center_0, in);
                vstorei_512(BV0[y - YSTART][0], v_center_0); 


                vloadi2_512(v_x_minus_0, BV0[y - YSTART + 2][0]);
                vloadi2_512(v_center_0,  BV1[y - YSTART + 2][0]);
                vloadi2_512(v_x_plus_0,  BV2[y - YSTART + 2][0]);
                Compute_1vector_512(v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_2_512(out, v_center_1, in);
                vstorei_512(BV0[y - YSTART + 1][0], v_center_1); 


                vloadi2_512(v_x_minus_1, BV0[y - YSTART + 3][0]);
                vloadi2_512(v_center_1,  BV1[y - YSTART + 3][0]);
                vloadi2_512(v_x_plus_1,  BV2[y - YSTART + 3][0]);
                Compute_1vector_512(v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_3_512(out, v_center_2, in);	
                vstorei_512(BV0[y - YSTART + 2][0], v_center_2); 

                vloadi2_512(v_x_minus_2, BV0[y - YSTART + 4][0]);
                vloadi2_512(v_center_2,  BV1[y - YSTART + 4][0]);
                vloadi2_512(v_x_plus_2,  BV2[y - YSTART + 4][0]);
                Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_4_512(out, v_center_0, in);
                vstorei_512(BV0[y - YSTART + 3][0], v_center_0); 

                vloadi2_512(v_x_minus_0, BV0[y - YSTART + 5][0]);
                vloadi2_512(v_center_0,  BV1[y - YSTART + 5][0]);
                vloadi2_512(v_x_plus_0,  BV2[y - YSTART + 5][0]);
                Compute_1vector_512(v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_5_512(out, v_center_1, in);
                vstorei_512(BV0[y - YSTART + 4][0], v_center_1); 

                vloadi2_512(v_x_minus_1, BV0[y - YSTART + 6][0]);
                vloadi2_512(v_center_1,  BV1[y - YSTART + 6][0]);
                vloadi2_512(v_x_plus_1,  BV2[y - YSTART + 6][0]);
                Compute_1vector_512(v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_6_512(out, v_center_2, in);	
                vstorei_512(BV0[y - YSTART + 5][0], v_center_2); 


                vloadi2_512(v_x_minus_2, BV0[y - YSTART + 7][0]);
                vloadi2_512(v_center_2,  BV1[y - YSTART + 7][0]);
                vloadi2_512(v_x_plus_2,  BV2[y - YSTART + 7][0]);
                Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_7_512(out, v_center_0, in);
                vstorei_512(BV0[y - YSTART + 6][0], v_center_0); 


                vloadi2_512(v_x_minus_0, BV0[y - YSTART + 8][0]);
                vloadi2_512(v_center_0,  BV1[y - YSTART + 8][0]);
                vloadi2_512(v_x_plus_0,  BV2[y - YSTART + 8][0]);
                Compute_1vector_512(v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_8_512(out, v_center_1, in);
                vstorei_512(BV0[y - YSTART + 7][0], v_center_1); 

                vloadi2_512(v_x_minus_1, BV0[y - YSTART + 9][0]);
                vloadi2_512(v_center_1,  BV1[y - YSTART + 9][0]);
                vloadi2_512(v_x_plus_1,  BV2[y - YSTART + 9][0]);
                Compute_1vector_512(v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_9_512(out, v_center_2, in);	
                vstorei_512(BV0[y - YSTART + 8][0], v_center_2); 

                vloadi2_512(v_x_minus_2, BV0[y - YSTART + 10][0]);
                vloadi2_512(v_center_2,  BV1[y - YSTART + 10][0]);
                vloadi2_512(v_x_plus_2,  BV2[y - YSTART + 10][0]);
                Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_10_512(out, v_center_0, in);
                vstorei_512(BV0[y - YSTART + 9][0], v_center_0); 

                vloadi2_512(v_x_minus_0, BV0[y - YSTART + 11][0]);
                vloadi2_512(v_center_0,  BV1[y - YSTART + 11][0]);
                vloadi2_512(v_x_plus_0,  BV2[y - YSTART + 11][0]);
                Compute_1vector_512(v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_11_512(out, v_center_1, in);
                vstorei_512(BV0[y - YSTART + 10][0], v_center_1); 

                vloadi2_512(v_x_minus_1, BV0[y - YSTART + 12][0]);
                vloadi2_512(v_center_1,  BV1[y - YSTART + 12][0]);
                vloadi2_512(v_x_plus_1,  BV2[y - YSTART + 12][0]);
                Compute_1vector_512(v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_12_512(out, v_center_2, in);	
                vstorei_512(BV0[y - YSTART + 11][0], v_center_2); 

                vloadi2_512(v_x_minus_2, BV0[y - YSTART + 13][0]);
                vloadi2_512(v_center_2,  BV1[y - YSTART + 13][0]);
                vloadi2_512(v_x_plus_2,  BV2[y - YSTART + 13][0]);
                Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_13_512(out, v_center_0, in);
                vstorei_512(BV0[y - YSTART + 12][0], v_center_0); 

                vloadi2_512(v_x_minus_0, BV0[y - YSTART + 14][0]);
                vloadi2_512(v_center_0,  BV1[y - YSTART + 14][0]);
                vloadi2_512(v_x_plus_0,  BV2[y - YSTART + 14][0]);
                Compute_1vector_512(v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_14_512(out, v_center_1, in);
                vstorei_512(BV0[y - YSTART + 13][0], v_center_1); 

                vloadi2_512(v_x_minus_1, BV0[y - YSTART + 15][0]);
                vloadi2_512(v_center_1,  BV1[y - YSTART + 15][0]);
                vloadi2_512(v_x_plus_1,  BV2[y - YSTART + 15][0]);
                Compute_1vector_512(v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0,
                                     v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_15_512(out, v_center_2, in);	
                vstorei_512(BV0[y - YSTART + 14][0], v_center_2); 



                if(y + VECLEN_INT_512 <= NY + YSTART - VECLEN_INT_512){
                    vloadi2_512(v_x_minus_2, BV0[y - YSTART + 16][0]);
                    vloadi2_512(v_center_2,  BV1[y - YSTART + 16][0]);
                    vloadi2_512(v_x_plus_2,  BV2[y - YSTART + 16][0]);
                } else {
                    vloadseti_512(v_x_minus_2, B, t, x-XSLOPE, y+VECLEN_INT_512);  
                    vloadseti_512(v_center_2, B, t, x, y+VECLEN_INT_512);  
                    vloadseti_512(v_x_plus_2, B, t, x+XSLOPE, y+VECLEN_INT_512); 
                }
                Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0); 
                Input_Output_i_16_512(out, v_center_0, in);
                vstorei_512(BV0[y - YSTART + 15][0], v_center_0);

                v_x_minus_0 = v_x_minus_1;
                v_center_0 = v_center_1;
                v_x_plus_0 = v_x_plus_1;
                v_x_minus_1 = v_x_minus_2;
                v_center_1 = v_center_2;
                v_x_plus_1 = v_x_plus_2;


                vstorei_512(B[(t)%2][x][y], out);

            }
         
            for ( y += 1; y <= NY+YSTART; y++){
                vloadseti_512(v_x_minus_2, B, t, x-XSLOPE, y);  
                vloadseti_512(v_center_2, B, t, x, y);  
                vloadseti_512(v_x_plus_2, B, t, x+XSLOPE, y); 

                Compute_1vector_512(v_x_minus_2, v_center_2, v_x_plus_2,
                                     v_x_minus_1, v_center_1, v_x_plus_1,
                                     v_x_minus_0, v_center_0, v_x_plus_0); 
                vstorei_512(tmp[0], v_center_0);


		    	B[(t+1)%2][x + STRIDE * 15][y - 1] = tmp[0];
		    	B[t%2    ][x + STRIDE * 14][y - 1] = tmp[1];
		    	B[(t+1)%2][x + STRIDE * 13][y - 1] = tmp[2];
		    	B[t%2    ][x + STRIDE * 12][y - 1] = tmp[3];
		    	B[(t+1)%2][x + STRIDE * 11][y - 1] = tmp[4];
		    	B[t%2    ][x + STRIDE * 10][y - 1] = tmp[5];
		    	B[(t+1)%2][x + STRIDE * 9 ][y - 1] = tmp[6];
		    	B[t%2    ][x + STRIDE * 8 ][y - 1] = tmp[7];
		    	B[(t+1)%2][x + STRIDE * 7 ][y - 1] = tmp[8];
		    	B[t%2    ][x + STRIDE * 6 ][y - 1] = tmp[9];
		    	B[(t+1)%2][x + STRIDE * 5 ][y - 1] = tmp[10];
		    	B[t%2    ][x + STRIDE * 4 ][y - 1] = tmp[11];
		    	B[(t+1)%2][x + STRIDE * 3 ][y - 1] = tmp[12];
		    	B[t%2    ][x + STRIDE * 2 ][y - 1] = tmp[13];
		    	B[(t+1)%2][x + STRIDE * 1 ][y - 1] = tmp[14];
		    	B[t%2    ][x             ][y - 1] = tmp[15];

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
			for ( y = YSTART; y <= NY + YSTART - VECLEN_INT_512; y += VECLEN_INT_512 ){

				vloadi2_512(v15, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 15][0]);
				vloadi2_512(v14, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 14][0]);
				vloadi2_512(v13, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 13][0]);
				vloadi2_512(v12, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 12][0]);
				vloadi2_512(v11, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 11][0]);
				vloadi2_512(v10, Btmp[x - (xx - XSLOPE)][(y - YSTART) + 10][0]);
				vloadi2_512(v9,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 9 ][0]);
				vloadi2_512(v8,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 8 ][0]);
				vloadi2_512(v7,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 7 ][0]);
				vloadi2_512(v6,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 6 ][0]);
				vloadi2_512(v5,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 5 ][0]);
				vloadi2_512(v4,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 4 ][0]);
				vloadi2_512(v3,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 3 ][0]);
				vloadi2_512(v2,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 2 ][0]);
				vloadi2_512(v1,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 1 ][0]);
				vloadi2_512(v0,  Btmp[x - (xx - XSLOPE)][(y - YSTART) + 0 ][0]);

                transposei_512(v0, v1, v2, v3, v4, v5, v6, v7,
                               v8, v9, v10, v11, v12, v13, v14, v15,
                               v0, v1, v2, v3, v4, v5, v6, v7,
                               v8, v9, v10, v11, v12, v13, v14, v15);

				vstorei_512( B[t%2    ][x + STRIDE * 15][y], v0);
				vstorei_512( B[(t+1)%2][x + STRIDE * 14][y], v1);
				vstorei_512( B[t%2    ][x + STRIDE * 13][y], v2);
				vstorei_512( B[(t+1)%2][x + STRIDE * 12][y], v3);
				vstorei_512( B[t%2    ][x + STRIDE * 11][y], v4);
				vstorei_512( B[(t+1)%2][x + STRIDE * 10][y], v5);
				vstorei_512( B[t%2    ][x + STRIDE * 9 ][y], v6);
				vstorei_512( B[(t+1)%2][x + STRIDE * 8 ][y], v7);
				vstorei_512( B[t%2    ][x + STRIDE * 7 ][y], v8);
				vstorei_512( B[(t+1)%2][x + STRIDE * 6 ][y], v9);
				vstorei_512( B[t%2    ][x + STRIDE * 5 ][y], v10);
				vstorei_512( B[(t+1)%2][x + STRIDE * 4 ][y], v11);
				vstorei_512( B[t%2    ][x + STRIDE * 3 ][y], v12);
				vstorei_512( B[(t+1)%2][x + STRIDE * 2 ][y], v13);
				vstorei_512( B[t%2    ][x + STRIDE * 1 ][y], v14);
				vstorei_512( B[(t+1)%2][x            ][y], v15);
			}
		} 
    

        for( t = tt ; t < tt + VECLEN_INT_512 ; t++){
			for ( x = xx + STRIDE * (VECLEN_INT_512 - 1 - (t - tt)); x < NX + XSTART; x++) {
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
    free_extra_array_512(AV);
}