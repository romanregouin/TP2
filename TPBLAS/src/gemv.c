#include "mnblas.h"

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,const MNCBLAS_TRANSPOSE TransA, const int M, const int N,const float alpha, const float *A, const int lda,
const float *X, const int incX, const float beta,float *Y, const int incY){
    register unsigned int i ;
    register unsigned int j ;
    float tmp;
    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && TransA != MNCblasNoTrans)) {
        for(i=0;i<M;i+=incY){
            for(j=0,tmp=0; (j < N); j+=incX ){
                tmp+=*(A+i*lda+j)*(*(X+j));
            }
            *(Y+i)=*(Y+i)*beta+tmp*alpha;
        }      
    }else{
        for(i=0;i<N;i++){
            for(j=0,tmp=0; (j < M); j+=incX  ){
                tmp+=*(A+i+j*lda)*(*(X+j));
            }
            *(Y+i)=*(Y+i)*beta+tmp*alpha;
        }
    }
    return;
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,MNCBLAS_TRANSPOSE TransA, const int M, const int N,const double alpha, const double *A, const int lda,
const double *X, const int incX, const double beta,double *Y, const int incY){
    register unsigned int i ;
    register unsigned int j ;
    double tmp;
    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && TransA != MNCblasNoTrans)) {
        for(i=0;i<M;i+=incY){
            for(j=0,tmp=0; (j < N); j+=incX ){
                tmp+=*(A+i*lda+j)*(*(X+j));
            }
            *(Y+i)=*(Y+i)*beta+tmp*alpha;
        }      
    }else{
        for(i=0;i<N;i++){
            for(j=0,tmp=0; (j < M); j+=incX ){
                tmp+=*(A+i+j*lda)*(*(X+j));
            }
            *(Y+i)=*(Y+i)*beta+tmp*alpha;
        }
    }
    return;
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,MNCBLAS_TRANSPOSE TransA, const int M, const int N, const void *alpha, const void *A, const int lda,
const void *X, const int incX, const void *beta,void *Y, const int incY){
    register unsigned int i ;
    register unsigned int j ;
    float im,re;
    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && TransA != MNCblasNoTrans)) {
        for(i=0;i<M;i+=incY){
            for(j=0,re=0,im=0; (j < N); j+=incX ){
                re+=*(((float*)A)+i*lda*2+j*2)*(*(((float*)X)+j*2)) - *(((float*)A)+i*lda*2+j*2+1)*(*(((float*)X)+j*2+1));
                im+=*(((float*)A)+i*lda*2+j*2)*(*(((float*)X)+j*2+1))+(*(((float*)X)+j*2))* (*(((float*)A)+i*lda*2+j*2+1));
            }
            *(((float*)Y)+i*2)=*(((float*)Y)+i*2)*(*((float*)beta))-*(((float*)Y)+i*2+1)*(*((float*)beta)+1)+ re*(*((float*)alpha))-im*(*((float*)alpha)+1);
            *(((float*)Y)+i*2+1)=*(((float*)Y)+i*2)*(*((float*)beta)+1)+*(((float*)Y)+i*2+1)*(*((float*)beta))+ re*(*((float*)alpha)+1)+im*(*((float*)alpha));
        }      
    }else{
        for(i=0;i<N;i++){
            for(j=0,re=0,im=0; (j < M); j+=incX ){
                re+=*(((float*)A)+i*2+j*2*lda)*(*(((float*)X)+j*2)) - *(((float*)A)+i*2+j*2*lda+1)*(*(((float*)X)+j*2+1));
                im+=*(((float*)A)+i*2+j*2*lda)*(*(((float*)X)+j*2+1))+(*(((float*)X)+j*2))* (*(((float*)A)+i*2+j*2*lda+1));
            }
            *(((float*)Y)+i*2)=*(((float*)Y)+i*2)*(*((float*)beta))-*(((float*)Y)+i*2+1)*(*((float*)beta)+1)+ re*(*((float*)alpha))-im*(*((float*)alpha)+1);
            *(((float*)Y)+i*2+1)=*(((float*)Y)+i*2)*(*((float*)beta)+1)+*(((float*)Y)+i*2+1)*(*((float*)beta))+ re*(*((float*)alpha)+1)+im*(*((float*)alpha));

        }
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,MNCBLAS_TRANSPOSE TransA, const int M, const int N,const void *alpha, const void *A, const int lda,
const void *X, const int incX, const void *beta,void *Y, const int incY){
    register unsigned int i ;
    register unsigned int j ;
    double im,re;
    if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) || (layout == MNCblasColMajor && TransA != MNCblasNoTrans)) {
        for(i=0;i<M;i+=incY){
            for(j=0,re=0,im=0; (j < N); j+=incX ){
                re+=*(((double*)A)+i*lda*2+j*2)*(*(((double*)X)+j*2)) - *(((double*)A)+i*lda*2+j*2+1)*(*(((double*)X)+j*2+1));
                im+=*(((double*)A)+i*lda*2+j*2)*(*(((double*)X)+j*2+1))+(*(((double*)X)+j*2))* (*(((double*)A)+i*lda*2+j*2+1));
            }
            *(((double*)Y)+i*2)=*(((double*)Y)+i*2)*(*((double*)beta))-*(((double*)Y)+i*2+1)*(*((double*)beta)+1)+ re*(*((double*)alpha))-im*(*((double*)alpha)+1);
            *(((double*)Y)+i*2+1)=*(((double*)Y)+i*2)*(*((double*)beta)+1)+*(((double*)Y)+i*2+1)*(*((double*)beta))+ re*(*((double*)alpha)+1)+im*(*((double*)alpha));
        }      
    }else{
        for(i=0;i<N;i++){
            for(j=0,re=0,im=0; (j < M); j+=incX ){
                re+=*(((double*)A)+i*2+j*2*lda)*(*(((double*)X)+j*2)) - *(((double*)A)+i*2+j*2*lda+1)*(*(((double*)X)+j*2+1));
                im+=*(((double*)A)+i*2+j*2*lda)*(*(((double*)X)+j*2+1))+(*(((double*)X)+j*2))* (*(((double*)A)+i*2+j*2*lda+1));
            }
            *(((double*)Y)+i*2)=*(((double*)Y)+i*2)*(*((double*)beta))-*(((double*)Y)+i*2+1)*(*((double*)beta)+1)+ re*(*((double*)alpha))-im*(*((double*)alpha)+1);
            *(((double*)Y)+i*2+1)=*(((double*)Y)+i*2)*(*((double*)beta)+1)+*(((double*)Y)+i*2+1)*(*((double*)beta))+ re*(*((double*)alpha)+1)+im*(*((double*)alpha));

        }
    }
}