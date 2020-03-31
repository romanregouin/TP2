#include "mnblas.h"
#include "stdlib.h"
#include "complexe2.h"

//Matrtice Carré / lda inutile + faire attention taille des matrice a ce quelle ne dépasse la taille du cache 256*256 ou 512*512
//vecteur 1000 4000 elements
//PAS DE TRANSPOSE

inline void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
     MNCBLAS_TRANSPOSE TransB, const int M, const int N,
     const int K, const float alpha, const float *A,
     const int lda, const float *B, const int ldb,
     const float beta, float *C, const int ldc){
    
    //int newrowA, newrowB, newcolA, newcolB;
    //applyOperation(TransA, layout, &newrowA, &newcolA, A, M, K); //op(A)
    //applyOperation(TransB, layout, &newrowB, &newcolB, B, K, N); //op(B)
    float* tmp = MultiplyByScalarConstMat(alpha,A,M,K);
    float* tmp1 = ProduitMatriciel(tmp,B,M,K,N);
    MultiplyByScalar(beta,C,M,N);
    AdditionMatriciel(tmp1,C,M,N);
    free(tmp);
    free(tmp1);
    return;
}

inline void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
     MNCBLAS_TRANSPOSE TransB, const int M, const int N,
     const int K, const double alpha, const double *A,
     const int lda, const double *B, const int ldb,
     const double beta, double *C, const int ldc){

    double* tmp = MultiplyByScalarConstMatDouble(alpha,A,M,K);
    double* tmp1 = ProduitMatricielDouble(tmp,B,M,K,N);
    MultiplyByScalarDouble(beta,C,M,N);
    AdditionMatricielDouble(tmp1,C,M,N);
    free(tmp);
    free(tmp1);
}

inline void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
     MNCBLAS_TRANSPOSE TransB, const int M, const int N,
     const int K, const void *alpha, const void *A,
     const int lda, const void *B, const int ldb,
     const void *beta, void *C, const int ldc){

    A = (complexe_float_t**)A;
    B = (complexe_float_t**)B;
    C = (complexe_float_t**)C;
    complexe_float_t alpha_complexe;
    alpha_complexe.real = ((float*)alpha)[0];
    alpha_complexe.imaginary = ((float*)alpha)[1]; 
    complexe_float_t beta_complexe;
    beta_complexe.real = ((float*)beta)[0];
    beta_complexe.imaginary = ((float*)beta)[1]; 
    void* tmp = MultiplyByScalarConstMatComplexe(alpha_complexe,A,M,K);
    void* tmp1 = ProduitMatricielComplexe(tmp,B,M,K,N);
    MultiplyByScalarComplexe(beta_complexe,C,M,N);
    AdditionMatricielComplexe(tmp1,C,M,N);
    free(tmp);
    free(tmp1);
}

inline void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
     MNCBLAS_TRANSPOSE TransB, const int M, const int N,
     const int K, const void *alpha, const void *A,
     const int lda, const void *B, const int ldb,
     const void *beta, void *C, const int ldc){
    
    A = (complexe_double_t**)A;
    B = (complexe_double_t**)B;
    C = (complexe_double_t**)C;
    complexe_double_t alpha_complexe;
    alpha_complexe.real = ((double*)alpha)[0];
    alpha_complexe.imaginary = ((double*)alpha)[1]; 
    complexe_double_t beta_complexe;
    beta_complexe.real = ((double*)beta)[0];
    beta_complexe.imaginary = ((double*)beta)[1]; 
    void* tmp = MultiplyByScalarConstMatComplexeDouble(alpha_complexe,A,M,K);
    void* tmp1 = ProduitMatricielComplexeDouble(tmp,B,M,K,N);
    MultiplyByScalarComplexeDouble(beta_complexe,C,M,N);
    AdditionMatricielComplexeDouble(tmp1,C,M,N);
    free(tmp);
    free(tmp1);
}

inline float* applyOperation(MNCBLAS_TRANSPOSE op, MNCBLAS_LAYOUT layout, int* newrow, int* newcol, const float* m, const int row, const int col) {
    float* tmp = malloc(col*row*sizeof(float));
    int pos = 0;
    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            tmp[pos] = m[i+j*col];
            pos++;
        }
    }
    *newcol = row;
    *newrow = col;
    return tmp;
}

inline float* MultiplyByScalarConstMat(const float alpha, const float* m, int row, int col){
    int lg = row*col;
    float* tmp = malloc(row*col*sizeof(float));
    for(int i=0;i<lg;i++){
        tmp[i] = m[i] * alpha; 
    }
    return tmp;
}

inline double* MultiplyByScalarConstMatDouble(const double alpha, const double* m, int row, int col){
    int lg = row*col;
    double* tmp = malloc(row*col*sizeof(double));
    for(int i=0;i<lg;i++){
        tmp[i] = m[i] * alpha; 
    }
    return tmp;
}

inline complexe_float_t* MultiplyByScalarConstMatComplexe(complexe_float_t alpha, const complexe_float_t* m, int row, int col){
    int lg = row*col;
    complexe_float_t tmp;
    complexe_float_t* tmp_res = malloc(row*col*sizeof(complexe_float_t));
    for(int i=0;i<lg;i++){
        tmp.imaginary = m[i].real * alpha.imaginary + m[i].imaginary * alpha.real;
        tmp.real = m[i].real * alpha.real - (m[i].imaginary * alpha.imaginary);
        tmp_res[i].real = tmp.real;
        tmp_res[i].imaginary = tmp.imaginary;
    }
    return tmp_res;
}

inline complexe_double_t* MultiplyByScalarConstMatComplexeDouble(complexe_double_t alpha, const complexe_double_t* m, int row, int col){
    int lg = row*col;
    complexe_double_t tmp;
    complexe_double_t* tmp_res = malloc(row*col*sizeof(complexe_double_t));
    for(int i=0;i<lg;i++){
        tmp.imaginary = m[i].real * alpha.imaginary + m[i].imaginary * alpha.real;
        tmp.real = m[i].real * alpha.real - (m[i].imaginary * alpha.imaginary);
        tmp_res[i].real = tmp.real;
        tmp_res[i].imaginary = tmp.imaginary;
    }
    return tmp_res;
}


inline void MultiplyByScalar(const float alpha, float* m, int row, int col){
    int lg = row*col;
    for(int i=0;i<lg;i++){
        m[i] *= alpha; 
    }
}

inline void MultiplyByScalarDouble(const double alpha, double* m, int row, int col){
    int lg = row*col;
    for(int i=0;i<lg;i++){
        m[i] *= alpha; 
    }
}

inline void MultiplyByScalarComplexe(complexe_float_t alpha, complexe_float_t* m, int row, int col){
    int lg = row*col;
    complexe_float_t tmp;
    for(int i=0;i<lg;i++){
        tmp.imaginary = m[i].real * alpha.imaginary + m[i].imaginary * alpha.real;
        tmp.real = m[i].real * alpha.real - (m[i].imaginary * alpha.imaginary);
        m[i].real = tmp.real;
        m[i].imaginary = tmp.imaginary;
    }
}

inline void MultiplyByScalarComplexeDouble(complexe_double_t alpha, complexe_double_t* m, int row, int col){
    int lg = row*col;
    complexe_double_t tmp;
    for(int i=0;i<lg;i++){
        tmp.imaginary = m[i].real * alpha.imaginary + m[i].imaginary * alpha.real;
        tmp.real = m[i].real * alpha.real - (m[i].imaginary * alpha.imaginary);
        m[i].real = tmp.real;
        m[i].imaginary = tmp.imaginary;
    }
}

inline float* ProduitMatriciel(const float* a, const float* b,const int m, const int k, const int n){
    float* res = malloc(m*n*sizeof(float));
    float result = 0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int z=0; z<k;z++){
                result += a[z+i*k]*b[j+z*n];
            }
            res[j+i*n] = result;
            result = 0;
        }
    }
    return res;
}

inline double* ProduitMatricielDouble(const double* a, const double* b,const int m, const int k, const int n){
    double* res = malloc(m*n*sizeof(double));
    double result = 0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int z=0; z<k;z++){
                result += a[z+i*k]*b[j+z*n];
            }
            res[j+i*n] = result;
            result = 0;
        }
    }
    return res;
}

inline complexe_float_t* ProduitMatricielComplexe(const complexe_float_t* a, const complexe_float_t* b,const int m, const int k, const int n){
    complexe_float_t* res = malloc(m*n*sizeof(complexe_float_t));
    complexe_float_t result;
    result.real = 0;
    result.imaginary = 0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int z=0; z<k;z++){
                result.imaginary += a[z+i*k].real * b[j+z*n].imaginary + a[z+i*k].imaginary * b[j+z*n].real;
                result.real += a[z+i*k].real * b[j+z*n].real - (a[z+i*k].imaginary * b[j+z*n].imaginary);
            }
            res[j+i*n] = result;
            result.real = 0;
            result.imaginary = 0;
        }
    }
    return res;
}

inline complexe_double_t* ProduitMatricielComplexeDouble(const complexe_double_t* a, const complexe_double_t* b,const int m, const int k, const int n){
    complexe_double_t* res = malloc(m*n*sizeof(complexe_double_t));
    complexe_double_t result;
    result.real = 0;
    result.imaginary = 0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int z=0; z<k;z++){
                result.imaginary += a[z+i*k].real * b[j+z*n].imaginary + a[z+i*k].imaginary * b[j+z*n].real;
                result.real += a[z+i*k].real * b[j+z*n].real - (a[z+i*k].imaginary * b[j+z*n].imaginary);
            }
            res[j+i*n] = result;
            result.real = 0;
            result.imaginary = 0;
        }
    }
    return res;
}

inline void AdditionMatriciel(float* a, float* b, const int m, const int n){
    for(int i=0;i<m*n;i++){
        b[i] = a[i] + b[i];
    }
}

inline void AdditionMatricielDouble(double* a, double* b, const int m, const int n){
    for(int i=0;i<m*n;i++){
        b[i] = a[i] + b[i];
    }
}

inline void AdditionMatricielComplexe(complexe_float_t* a, complexe_float_t* b, const int m, const int n){
    for(int i=0;i<m*n;i++){
        b[i].real = a[i].real + b[i].real;
        b[i].imaginary = a[i].imaginary + b[i].imaginary;
    }
}

inline void AdditionMatricielComplexeDouble(complexe_double_t* a, complexe_double_t* b, const int m, const int n){
    for(int i=0;i<m*n;i++){
        b[i].real = a[i].real + b[i].real;
        b[i].imaginary = a[i].imaginary + b[i].imaginary;
    }
}