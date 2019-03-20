#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct HeatTransfer heatTransfer;

struct HeatTransfer{
	int n;
	float BFL;
	float final_temp;
	float*** input_temp;
	float*** output_temp;
};

void printTemp(struct HeatTransfer ht){
	
	for (int i = 0; i < ht.n+2; ++i){
                for (int j = 0; j < ht.n+2; ++j){
                        for (int k = 0; k < ht.n+2; ++k){
				fprintf(stderr, " %f ", ht.input_temp[i][j][k]);
                        }
                        fprintf(stderr, "\n");
                }
        }
}

void Advance(struct HeatTransfer ht){
	for (int i = 1; i < ht.n+1; ++i){
		for (int j = 1; j < ht.n+1; j++){
			for (int k = 1; k < ht.n+1; k++){
				ht.output_temp[i][j][k] = (ht.input_temp[i-1][j][k] + ht.input_temp[i+1][j][k] + ht.input_temp[i][j-1][k] + ht.input_temp[i][j+1][k] + ht.input_temp[i][j][k-1] + ht.input_temp[i][j][k+1]) / 6.0;
			}
		}
	}
	for (int i = 1; i < ht.n+1; ++i){
		for (int j = 1; j < ht.n+1; j++){
			for (int k = 1; k < ht.n+1; k++){
				ht.input_temp[i][j][k] = ht.output_temp[i][j][k];
			}
		}	
	}
}

float calculateFinal_Temp(struct HeatTransfer ht){
	for (int i = 1; i < ht.n+1; ++i){
		for (int j = 1; j < ht.n+1; j++){
			for (int k = 1; k < ht.n+1; k++){
				ht.final_temp += ht.input_temp[i][j][k];			
			}
		}
	}
	fprintf(stderr, "Final Temp: %f\n", ht.final_temp);
}


float lerp(float x, float x1, float x2, float q00, float q01){
	if (x2 == x1){
		return q00;
	} else {
		return ((x2 - x) / (x2 - x1)) * q00 + ((x - x1) / (x2 - x1)) * q01;
	}
}

float getHeat(struct HeatTransfer ht, float x, float y, float z){
	float t1, t2;
	float xt = floor(x);	
	float yt = floor(y);	
	float zt = floor(z);	
	t1 = ht.input_temp[(int)xt][(int)yt][(int)zt];
	t2 = ht.input_temp[(int)ceil(x)][(int)floor(y)][(int)floor(z)];
	float x00 = lerp(x, floor(x), ceil(x), t1, t2);
	
	t1 = ht.input_temp[(int)floor(x)][(int)ceil(y)][(int)floor(z)];
	t2 = ht.input_temp[(int)ceil(x)][(int)ceil(y)][(int)floor(z)];
	float x01 = lerp(x, floor(x), ceil(x), t1, t2);
	
	float y0 = lerp(y, floor(y), ceil(y), x00, x01);
	
	t1 = ht.input_temp[(int)floor(x)][(int)floor(y)][(int)ceil(z)];
	t2 = ht.input_temp[(int)ceil(x)][(int)floor(y)][(int)ceil(z)];
	float x10 = lerp(x, floor(x), ceil(x), t1, t2);

	t1 = ht.input_temp[(int)floor(x)][(int)ceil(y)][(int)ceil(z)];
	t2 = ht.input_temp[(int)ceil(x)][(int)ceil(y)][(int)ceil(z)];
	float x11 = lerp(x, floor(x), ceil(x), t1, t2);

	float y1 = lerp(y, floor(y), ceil(y), x10, x11);

	float z0 = lerp(z, floor(z), ceil(z), y0, y1);
	
	fprintf(stderr, "z = %f\n", z0);
	return z0;
}


int main(void)
{
	int n_iters = 100;
	struct HeatTransfer ht;

	ht.n = 10;
	ht.input_temp = (float***) malloc((ht.n+2) * sizeof(float**));	
	ht.output_temp = (float***) malloc((ht.n+2) * sizeof(float**));	
	for (int i = 0; i < ht.n+2; i++){
		ht.input_temp[i] = (float**) malloc((ht.n+2) * sizeof(float*));
		ht.output_temp[i] = (float**) malloc((ht.n+2) * sizeof(float*));
		for (int j = 0; j < ht.n+2; j++){
			ht.input_temp[i][j] = (float*) malloc((ht.n+2) * sizeof(float));
			ht.output_temp[i][j] = (float*) malloc((ht.n+2) * sizeof(float));
		}
	}

	for (int i = 0; i < ht.n+2; ++i){
		for (int j = 0; j < ht.n+2; ++j){
			for (int k = 0; k < ht.n+2; ++k){
				if (i == 0 || j == 0 || k == 0 || i == ht.n+1 || j == ht.n+1 || k == ht.n+1){
					ht.input_temp[i][j][k] = 100.0;
					ht.output_temp[i][j][k] = 100.0;
				} else {
					ht.input_temp[i][j][k] = -100.0;
				}
			}
		}
	}
	clock_t t;
	t = clock();
	for (int i = 0; i < n_iters; i++){
		Advance(ht);
	}
	t = clock() - t;
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

	printf("Execution Time: %f seconds\n", time_taken);

	//printTemp(ht);
	//getHeat(ht, 1.0, 1.1, 1.2);
	calculateFinal_Temp(ht);

	for (int i = 0; i < ht.n+2; ++i){
		for (int j = 0; j < ht.n+2; ++j){
			free(ht.input_temp[i][j]);
			free(ht.output_temp[i][j]);
		}
		free(ht.input_temp[i]);
		free(ht.output_temp[i]);
	}
	free(ht.input_temp);
	free(ht.output_temp);
}
