#include "mnblas.h"

void mncblas_sswap(const int N, float *X, const int incX, float *Y, const int incY){

  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX, double *Y, const int incY){
  int i = 0;
  int j = 0;
  double save;
  for(;((i<N)&&(j<N));i+=incX,j+=incY){
    save = Y[j];
    Y[j] = X[i];
    X[i] = save;
  }
  return ;
}

void mncblas_cswap(const int N, void *X, const int incX, void *Y, const int incY){
  int i = 0;
  int j = 0;
  float save_rel;
  float save_img;
  for(;(i<N)&&(j<N);i+=incX,j+=incY){
    save_rel = ((float*)Y)[2*j];
    save_img = ((float*)Y)[2*j+1];
    ((float*)Y)[2*j] = ((float*)X)[2*j];
    ((float*)Y)[2*j+1] = ((float*)X)[2*j+1];
    ((float*)X)[2*j] = save_rel;
    ((float*)X)[2*j+1] = save_img;
  }
  return ;
}
void mncblas_zswap(const int N, void *X, const int incX, void *Y, const int incY){
  int i = 0;
  int j = 0;
  double save_rel;
  double save_img;
  for(;(i<N)&&(j<N);i+=incX,j+=incY){
    save_rel = ((double*)Y)[2*j];
    save_img = ((double*)Y)[2*j+1];
    ((double*)Y)[2*j] = ((double*)X)[2*j];
    ((double*)Y)[2*j+1] = ((double*)X)[2*j+1];
    ((double*)X)[2*j] = save_rel;
    ((double*)X)[2*j+1] = save_img;
  }
  return ;
}

