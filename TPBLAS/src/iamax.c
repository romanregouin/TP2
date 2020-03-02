#include "mnblas.h"

CBLAS_INDEX mnblas_isamax(const int N, const float  *X, const int incX){
    CBLAS_INDEX indice = 0;
    float max = X[indice];
    for(long unsigned int i=0;i<N;i+=incX){
        if(X[i]>max){
            max = X[i];
            indice = i;
        }
    }
    return indice;
}

CBLAS_INDEX mnblas_idamax(const int N, const double *X, const int incX){
    CBLAS_INDEX indice = 0;
    double max = X[indice];
    for(long unsigned int i=0;i<N;i+=incX){
        if(X[i]>max){
            max = X[i];
            indice = i;
        }
    }
    return indice;
}

CBLAS_INDEX mnblas_icamax(const int N, const void   *X, const int incX){
    CBLAS_INDEX indice = 0;
    float max = ((float*)X)[0] + ((float*)X)[1];
    for(long unsigned int i=0;i<N;i+=incX){
        float sum = ((float*)X)[2*i] + ((float*)X)[2*i+1];
        if(sum>max){
            max = sum;
            indice = i;
        }
    }
    return indice;
}

CBLAS_INDEX mnblas_izamax(const int N, const void   *X, const int incX){
    CBLAS_INDEX indice = 0;
    double max = ((double*)X)[0] + ((double*)X)[1];
    for(long unsigned int i=0;i<N;i+=incX){
        double sum = ((double*)X)[2*i] + ((double*)X)[2*i+1];
        if(sum>max){
            max = sum;
            indice = i;
        }
    }
    return indice;
}