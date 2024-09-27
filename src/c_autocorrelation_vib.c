#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

double* compute_autocorrelation(double *pos, int nostep, int maxtau){//Compute the autocorrelation of the tab pos

	double* res;
	double normalization = nostep;

	res = calloc(maxtau,sizeof(double));
	double number;

	for(int i=0 ; i<maxtau ; i++){
		number = 0;
		for(int j=i ; j<nostep ; j++){
			number = number + pos[j]*pos[j-i];
		}
		res[i]=number/normalization;
		normalization--;
	}
	return res;
}