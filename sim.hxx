#ifndef SIM_H
#define SIM_H

#include <iostream>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>

float lerp(float x, float x1, float x2, float q00, float q01){
    if (x2 == x1){
        return q00;
    } else {
        return ((x2 - x) / (x2 - x1)) * q00 + ((x - x1) / (x2 - x1)) * q01;
    }
}

class HeatTransfer
{
    public:
    
    
    int n;
    float coef;
    float ***input_temp;
    float ***output_temp;
    float final_temp;
    float BFL;
    
    public:
    HeatTransfer (int _n, float _BFL) {
        n = _n;
        BFL = _BFL;
        
        float c = 0.1;
        float d_t = 1/(4*c);
        
        coef = c * d_t;
        
        input_temp = (float ***) malloc(sizeof(float **) * (n+2));
        output_temp = (float ***) malloc(sizeof(float **) * (n+2));
        //#pragma omp parallel for
        for (int i = 0; i < n+2; i++) {
            input_temp[i] = (float **) malloc(sizeof(float *) * (n+2));
            output_temp[i] = (float **) malloc(sizeof(float *) * (n+2));
            for (int j = 0; j < n+2; j++) {
                input_temp[i][j] = (float *) malloc(sizeof(float) * (n+2));
                output_temp[i][j] = (float *) malloc(sizeof(float) * (n+2));
                for (int k = 0; k < n+2; k++) {
                    if (i == 0 || i == n+1 ||
                        j == 0 || j == n+1 ||
                        k == 0 || k == n+1) {
                        input_temp[i][j][k] = -100.0;
                        output_temp[i][j][k] = -100.0;
                    }
                    else {
                        input_temp[i][j][k] = 100.0;
                        output_temp[i][j][k] = 100.0;
                    }
                }
            }
        }
    }
    void Advance() {
        std::cerr << "Advancing" << std::endl;
        //for (int m = 0 ; m < 1; m++) {
        for (int i = 1; i < n+1; i++) {
            for (int j = 1; j < n+1; j++) {
                for (int k = 1; k < n+1; k++) {
                    output_temp[i][j][k] = input_temp[i][j][k] + coef * (input_temp[i-1][j][k] + input_temp[i+1][j][k] + input_temp[i][j-1][k] + input_temp[i][j+1][k] + input_temp[i][j][k-1] + input_temp[i][j][k+1] - 6.0 * input_temp[i][j][k]);
//                    std::cerr << output_temp[i][j][k] << std::endl;
                }
            }
        }
        for (int i = 1; i < n+1; i++) {
            for (int j = 1; j < n+1; j++) {
                for (int k = 1; k < n+1; k++) {
                    input_temp[i][j][k] = output_temp[i][j][k];
                }
            }
        }
        std::cerr << "Advanced" << std::endl;
        //}
    }
    
    float getHeat(float x, float y, float z) {
        float t1, t2;
        float xt = floor(x);
        float yt = floor(y);
        float zt = floor(z);
        t1 = input_temp[(int)xt][(int)yt][(int)zt];
        t2 = input_temp[(int)ceil(x)][(int)floor(y)][(int)floor(z)];
        float x00 = lerp(x, floor(x), ceil(x), t1, t2);
        
        t1 = input_temp[(int)floor(x)][(int)ceil(y)][(int)floor(z)];
        t2 = input_temp[(int)ceil(x)][(int)ceil(y)][(int)floor(z)];
        float x01 = lerp(x, floor(x), ceil(x), t1, t2);
        
        float y0 = lerp(y, floor(y), ceil(y), x00, x01);
        
        t1 = input_temp[(int)floor(x)][(int)floor(y)][(int)ceil(z)];
        t2 = input_temp[(int)ceil(x)][(int)floor(y)][(int)ceil(z)];
        float x10 = lerp(x, floor(x), ceil(x), t1, t2);
        
        t1 = input_temp[(int)floor(x)][(int)ceil(y)][(int)ceil(z)];
        t2 = input_temp[(int)ceil(x)][(int)ceil(y)][(int)ceil(z)];
        float x11 = lerp(x, floor(x), ceil(x), t1, t2);
        
        float y1 = lerp(y, floor(y), ceil(y), x10, x11);
        
        float z0 = lerp(z, floor(z), ceil(z), y0, y1);
        
        //fprintf(stderr, "z = %f\n", z0);
        //std::cerr << "Heating" << std::endl;
        return z0;
    }
};


#endif
