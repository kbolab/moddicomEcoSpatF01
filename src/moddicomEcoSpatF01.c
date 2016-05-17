#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void internalLoop(double *D, int *dimD, double *dStar, double *W){
  int N;
  N = *dimD;
  int i;
  double sum;
  
  for (i=0 ; i<N*N; i++){ 	
    if (D[i] >= *dStar || D[i] < 10e-7) { W[i] = 0.00;}
    else if (D[i] < *dStar && D[i] > 10e-7){ W[i] = pow(D[i],-2);}
  }
}