#include "mnblas.h"

float  mnblas_sasum(const int N, const float *X, const int incX){
    float res = 0;
    float val;
    for(int i=0;i<N;i+=incX){
        val = X[i];
        if(val<0){
            val = -val;
        }
        res += val;
    }
    return res;
}

double mnblas_dasum(const int N, const double *X, const int incX){
    double res = 0;
    double val;
    for(int i=0;i<N;i+=incX){
        val = X[i];
        if(val<0){
            val = -val;
        }
        res += val;
    }
    return res;
}

float  mnblas_scasum(const int N, const void *X, const int incX){
    float res = 0;
    float reel;
    float img;
    for(int i=0;i<N;i+=incX){
        reel = ((float*)X)[2*i];
        img = ((float*)X)[2*i+1];
        if(reel<0){
            reel = -reel;
        }
        if(img<0){
            img = -img;
        }
        res += reel + img;
    }
    return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX){
    double res = 0;
    double reel;
    double img;
    for(int i=0;i<N;i+=incX){
        reel = ((double*)X)[2*i];
        img = ((double*)X)[2*i+1];
        if(reel<0){
            reel = -reel;
        }
        if(img<0){
            img = -img;
        }
        res += reel + img;
    }
    return res;
}