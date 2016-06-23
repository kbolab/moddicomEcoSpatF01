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



void varMoranNull(double *Wstar, double *t_Wstar, int *dimWstar, double *S1, double *S2){

  int N;
  N = *dimWstar;
 //per S1 prendo la trasposta e sommo stesso indice
  int i;
  for (i=0 ; i<N*N; i++){
    *S1 = (*S1 + pow(Wstar[i] + t_Wstar[i],2));
     //printf("\n i=%d, Wij=%lf, Wji=%lf, S1=%lf ", i , Wstar[i], t_Wstar[i], *S1);
  }
  *S1 = *S1/2;
  //printf("\n  S1=%lf",*S1);


  // per S2 mi muovo sui primi N elementi di matrice e sua traposta --> poi il quadrato della somma --> poi la somma dei quadrati
  int j;
  double tmp;
  for (i=0 ; i<N; i++){
    for (j=0 ; j<N; j++){
    tmp = tmp + Wstar[i * N + j] + t_Wstar[i * N + j];
    //printf("\n i=%d, j=%d, W=%lf, S2=%lf, i * N + j=%d", i ,j, Wstar[i * N + j],*S2,i * N + j);
    }
    *S2 = *S2 + pow(tmp,2);
    tmp=0;
  }
  //printf("\n  S2=%lf",*S2);

 //Per S3 mi servono i valorio degli scarti dalla media (grayVar)

//  double tmp_gray4;
//  double tmp_gray2;
//  for (i=0 ; i<N*N; i++){
//    tmp_gray4 = (tmp_gray4 + pow(grayVar[i],4));
//    tmp_gray2 = (tmp_gray2 + pow(grayVar[i],2));
//    printf("\n tmp_gray4=%lf, tmp_gray2=%lf", tmp_gray4,tmp_gray2);
//  }
//  *S3 = N * (tmp_gray4) * (1/(pow(tmp_gray2,2))) ;
//  printf("\n S3=%lf", *S3);


  // qua mi devo solo calcolare il quadrato di ogni elemento di matrice e sommarlo --> poi calcolo S$ e S% con le formule
//  double tmpW2;
//  for (i=0 ; i<N*N; i++){
//    tmpW2 = tmpW2 + Wstar[i];
//    printf("\n i=%d, Wij=%lf, tmpW2=%lf ", i , Wstar[i], tmpW2);
//  }
//  *S4= ( pow(N,2) - 3*N + 3 ) * (*S1) - N * (*S2) +  3 * tmpW2;
//  *S5= ( pow(N,2) - N ) * (*S1) - 2 * N * (*S2) +  6 * tmpW2;
//  printf("\n S4=%lf, S5=%lf ", *S4, *S5);

 // Mi calcolo finalmente la varianza attesa per Hp nulla
 //  *varI = ( ( N * (*S4) - (*S3) * (*S5) ) / (( N - 1) * ( N - 2) * ( N - 3) * tmpW2) ) - pow(*expectI,2);
 //  printf("\n  varI=%lf",*varI);
}
