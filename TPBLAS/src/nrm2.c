#include "mnblas.h"
#include <math.h>

float  mnblas_snrm2(const int N, const float *X, const int incX){
    float res = 0;
    for(int i=0;i<N;i+=incX){
        res += X[i]*X[i];
    }
    return sqrtf(res);
}

double mnblas_dnrm2(const int N, const double *X, const int incX){
    double res = 0;
    for(int i=0;i<N;i+=incX){
        res += X[i]*X[i];
    }
    return sqrt(res);
}

float  mnblas_scnrm2(const int N, const void *X, const int incX){
    float somme_img = 0;
    float somme_real = 0; 
    for(int i=0;i<N;i+=incX){
        float real = ((float*)X)[2*i];
        float img = ((float*)X)[2*i+1];
        somme_img += img;
        somme_real += real;
    }
    return sqrtf(somme_real*somme_real+somme_img*somme_img);
}

double mnblas_dznrm2(const int N, const void *X, const int incX){
    double somme_img = 0;
    double somme_real = 0;
    for(int i=0;i<N;i+=incX){
        double real = ((double*)X)[2*i];
        double img = ((double*)X)[2*i+1];
        somme_img += img;
        somme_real += real;
    }
    return sqrt(somme_real*somme_real+somme_img*somme_img);
}