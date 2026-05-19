#include "define.h"


void naive_scalar(double* A, int NX, int NY, int NZ, int T) {
	double (* A_correct)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART] =  (double (*)[NX + 2 * XSTART][ NY + 2 * YSTART][ NZ + 2 * ZSTART]) A;
    int t, x, y, z;
    for (t = 0; t < T; t++) {
		for (x = XSTART; x < NX + XSTART; x++) {
			for (y = YSTART; y < NY + YSTART; y++) {
				#pragma novector
				for (z = ZSTART; z < NZ + ZSTART; z++) {
					Compute_scalar(A_correct,t,x,y,z);
				}
			}
		}
		//         if(t == 3 ) \
        //             printf("Vec: %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",   \
        //                    A_correct[(t+1)%2][8][1][4],\
        // /* mid */          A_correct[(t)%2  ][8][1][4], \
        // /* left */         A_correct[(t)%2  ][7][1][4], \
        // /* right */        A_correct[(t)%2  ][9][1][4], \
        // /* forward */      A_correct[(t)%2  ][8][1][5], \
        // /* backward */     A_correct[(t)%2  ][8][1][3], \
        // /* down */         A_correct[(t)%2  ][8][0][4], \
        // /* up */           A_correct[(t)%2  ][8][2][4]);  
            // if(t == 3 ) {
            //             x=5; y=1; z=4;
            //             printf("Vec: (%d,%d,%d) %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",x,y,z,   \
            //                    A_correct[(t+1)%2][x  ][y  ][z   ],\
            // /* mid */          A_correct[(t)%2  ][x  ][y  ][z   ], \
            // /* left */         A_correct[(t)%2  ][x-1][y  ][z   ], \
            // /* right */        A_correct[(t)%2  ][x+1][y  ][z   ], \
            // /* forward */      A_correct[(t)%2  ][x  ][y  ][z-1 ], \
            // /* backward */     A_correct[(t)%2  ][x  ][y  ][z+1 ], \
            // /* down */         A_correct[(t)%2  ][x  ][y-1][z   ], \
            // /* up */           A_correct[(t)%2  ][x  ][y+1][z   ]);  }

	}

}