#include "define.h"

void vectime(int* A, int NX, int NY, int T) {

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
        //vloadi(a)     返回一个 __m256i，外面自己赋值
        //vloadi2(b,a)  直接把 load 结果赋给 b
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
				vstorei( B[t%2    ][x + STRIDE *7][y], vthree);
				vstorei( B[(t+1)%2][x + STRIDE *6][y], v_x_minus_2);
				vstorei( B[t%2    ][x + STRIDE *5][y], v_x_minus_1);
				vstorei( B[(t+1)%2][x + STRIDE *4][y], v_x_minus_0);
				vstorei( B[t%2    ][x + STRIDE *3][y], vone);
				vstorei( B[(t+1)%2][x + STRIDE *2][y], vzero);
				vstorei( B[t%2    ][x + STRIDE *1][y], out);
				vstorei( B[(t+1)%2][x            ][y], in);
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


            vloadi2(v_x_minus_1, B[(t+1)%2  ][x-XSLOPE][y]);
            vloadi2(v_center_1, B[(t+1)%2  ][x       ][y]);
            vloadi2(v_x_plus_1, B[(t+1)%2  ][x+XSLOPE][y]);

            for ( y = YSTART; y <= NY + YSTART - VECLEN_INT; y+= VECLEN_INT) {

 
                vloadi2(in, B[(t)%2][x+STRIDE*VECLEN_INT][y]); 


                vloadi2(v_x_minus_2, B[(t)%2  ][x-XSLOPE+STRIDE  ][y]);
                vloadi2(v_center_2, B[(t)%2  ][x+STRIDE  ][y]); 
                vloadi2(v_x_plus_2, B[(t)%2  ][x+XSLOPE+STRIDE  ][y]); 
                Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_1(out, v_center_0, in);
                vstorei(B[(t+1)%2  ][x+STRIDE][y], v_center_0); 



                vloadi2(v_x_minus_0, B[(t+1)%2][x-XSLOPE+STRIDE*2][y]);
                vloadi2(v_center_0, B[(t+1)%2][x+STRIDE*2][y]);
                vloadi2(v_x_plus_0, B[(t+1)%2][x+XSLOPE+STRIDE*2][y]);
                Compute_1vector(v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_2(out, v_center_1, in);
                vstorei(B[(t)%2][x+STRIDE*2][y], v_center_1);


                vloadi2(v_x_minus_1, B[(t)%2][x-XSLOPE+STRIDE*3][y]);
                vloadi2(v_center_1, B[(t)%2][x       +STRIDE*3][y]);
                vloadi2(v_x_plus_1, B[(t)%2][x+XSLOPE+STRIDE*3][y]);
                Compute_1vector(v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_3(out, v_center_2, in);	
                vstorei(B[(t+1)%2  ][x+STRIDE*3][y], v_center_2);


                vloadi2(v_x_minus_2, B[(t+1)%2  ][x-XSLOPE+STRIDE*4  ][y]);
                vloadi2(v_center_2, B[(t+1)%2  ][x+STRIDE*4  ][y]); 
                vloadi2(v_x_plus_2, B[(t+1)%2  ][x+XSLOPE+STRIDE*4  ][y]); 
                Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_4(out, v_center_0, in);
                vstorei(B[(t)%2  ][x+STRIDE*4][y], v_center_0); 


                vloadi2(v_x_minus_0, B[(t)%2][x-XSLOPE+STRIDE*5][y]);
                vloadi2(v_center_0, B[(t)%2][x+STRIDE*5][y]);
                vloadi2(v_x_plus_0, B[(t)%2][x+XSLOPE+STRIDE*5][y]);
                Compute_1vector(v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_5(out, v_center_1, in);
                vstorei(B[(t+1)%2][x+STRIDE*5][y], v_center_1);


                vloadi2(v_x_minus_1, B[(t+1)%2][x-XSLOPE+STRIDE*6][y]);
                vloadi2(v_center_1, B[(t+1)%2][x       +STRIDE*6][y]);
                vloadi2(v_x_plus_1, B[(t+1)%2][x+XSLOPE+STRIDE*6][y]);
                Compute_1vector(v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2); 
                Input_Output_i_6(out, v_center_2, in);	
                vstorei(B[(t)%2  ][x+STRIDE*6][y], v_center_2);



                vloadi2(v_x_minus_2, B[(t)%2  ][x-XSLOPE+STRIDE*7  ][y]);
                vloadi2(v_center_2, B[(t)%2  ][x+STRIDE*7  ][y]); 
                vloadi2(v_x_plus_2, B[(t)%2  ][x+XSLOPE+STRIDE*7  ][y]); 
                Compute_1vector(v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1,\
                                v_x_minus_0, v_center_0, v_x_plus_0);
                Input_Output_i_7(out, v_center_0, in);
                vstorei(B[(t+1)%2  ][x+STRIDE*7][y], v_center_0); 



                if(y + VECLEN_INT <= NY + YSTART - VECLEN_INT){
                    vloadi2(v_x_minus_0, B[(t+1)%2][x-XSLOPE][y+VECLEN_INT]);
                    vloadi2(v_center_0, B[(t+1)%2][x       ][y+VECLEN_INT]);
                    vloadi2(v_x_plus_0, B[(t+1)%2][x+XSLOPE][y+VECLEN_INT]);
                } else {
                    vloadseti(v_x_minus_0, B, t, x-XSLOPE, y+VECLEN_INT);  
                    vloadseti(v_center_0, B, t, x, y+VECLEN_INT);  
                    vloadseti(v_x_plus_0, B, t, x+XSLOPE, y+VECLEN_INT); 
                }
                Compute_1vector(v_x_minus_0, v_center_0, v_x_plus_0,\
                                v_x_minus_2, v_center_2, v_x_plus_2,\
                                v_x_minus_1, v_center_1, v_x_plus_1); 
                Input_Output_i_8(out, v_center_1, in);
                vstorei(B[(t)%2][x+STRIDE*8][y], v_center_1);

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
        }

        for ( y = YSTART; y < NY + YSTART ; y++ ){
            B[1][XSTART - XSLOPE][y] = B[0][XSTART - XSLOPE][y];
        }

        xx = x;
        for(x = xx - XSLOPE; x <= xx +XSLOPE; x++){
			for ( y = YSTART; y <= NY + YSTART - VECLEN_INT; y += VECLEN_INT ){
				vloadi2(vtwo, B[t%2    ][x + STRIDE * 7][y]);
				vloadi2(v_x_plus_2, B[(t+1)%2][x + STRIDE * 6][y]);
				vloadi2(v_x_plus_1, B[t%2    ][x + STRIDE * 5][y]);
				vloadi2(v_x_plus_0, B[(t+1)%2][x + STRIDE * 4][y]);
				vloadi2(vrotatei_high2low, B[t%2    ][x + STRIDE * 3][y]);
				vloadi2(v_center_2, B[(t+1)%2][x + STRIDE * 2][y]);
				vloadi2(v_center_1, B[t%2    ][x + STRIDE * 1][y]);
				vloadi2(v_center_0, B[(t+1)%2][x             ][y]);
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
}
