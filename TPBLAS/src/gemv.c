#include "mnblas.h";

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,const MNCBLAS_TRANSPOSE TransA, const int M, const int N,const float alpha, const float *A, const int lda,
const float *X, const int incX, const float beta,float *Y, const int incY){
    if(TransA==MNCblasNoTrans){
        //y := alpha*A*x + beta*y;
        
    }else if(TransA==MNCblasTrans){
        //y := alpha*A'*x + beta*y;
    }else if(TransA==MNCblasConjTrans){
        //y := alpha *conjg(A')*x + beta*y.
    }
    return;
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,MNCBLAS_TRANSPOSE TransA, const int M, const int N,const double alpha, const double *A, const int lda,
const double *X, const int incX, const double beta,double *Y, const int incY){
    return;
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha, const void *A, const int lda,
const void *X, const int incX, const void *beta,void *Y, const int incY){
    return;
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,MNCBLAS_TRANSPOSE TransA, const int M, const int N,const void *alpha, const void *A, const int lda,
const void *X, const int incX, const void *beta,void *Y, const int incY){
    return;
}