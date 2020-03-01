#include "mnblas.h"

CBLAS_INDEX mnblas_isamin(const int N, const float  *X, const int incX){
    CBLAS_INDEX indice = 0;
    float min = X[indice];
    for(long unsigned int i=0;i<N;i+=incX){
        if(X[i]<min){
            min = X[i];
            indice = i;
        }
    }
    return indice;
}

CBLAS_INDEX mnblas_idamin(const int N, const double *X, const int incX){
    CBLAS_INDEX indice = 0;
    double min = X[indice];
    for(long unsigned int i=0;i<N;i+=incX){
        if(X[i]<min){
            min = X[i];
            indice = i;
        }
    }
    return indice;
}

CBLAS_INDEX mnblas_icamin(const int N, const void   *X, const int incX){
    CBLAS_INDEX indice = 0;
    float min = ((float*)X)[0] + ((float*)X)[1];
    for(long unsigned int i=0;i<N;i+=incX){
        float sum = ((float*)X)[2*i] + ((float*)X)[2*i+1];
        if(sum<min){
            min = sum;
            indice = i;
        }
    }
    return indice;
}

CBLAS_INDEX mnblas_izamin(const int N, const void   *X, const int incX){
    CBLAS_INDEX indice = 0;
    double min = ((double*)X)[0] + ((double*)X)[1];
    for(long unsigned int i=0;i<N;i+=incX){
        double sum = ((double*)X)[2*i] + ((double*)X)[2*i+1];
        if(sum<min){
            min = sum;
            indice = i;
        }
    }
    return indice;
}