#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#include "flop.h"

//#define VECSIZE    36//65536

#define NB_FOIS    10

float* matrice_float(int m, int n){
    float* res = malloc(m*n*sizeof(float));
    float z = 1.0;
    for(int i=0;i<m*n;i++){
        res[i] = z;
        z++;
    }
    return res;
}

void print_matrice_float(float* a, int m, int n){
    printf("Matrice :\n");
    for(int i=0;i<m*n;i++){
        if((i%n==0) && (i!=0)){
            printf("\n");
        }
        printf("%f  ",a[i]);
    }
    printf("\n");
}

double* matrice_double(int m, int n){
    double* res = malloc(m*n*sizeof(double));
    double z = 1.0;
    for(int i=0;i<m*n;i++){
        res[i] = z;
        z++;
    }
    return res;
}

void print_matrice_double(double* a, int m, int n){
    printf("Matrice :\n");
    for(int i=0;i<m*n;i++){
        if((i%n==0) && (i!=0)){
            printf("\n");
        }
        printf("%lf  ",a[i]);
    }
    printf("\n");
}

complexe_float_t* matrice_complexe_float_t(int m, int n){
    complexe_float_t* res = malloc(m*n*sizeof(complexe_float_t));
    complexe_float_t z;
    z.real = 1.0;
    z.imaginary = 1.0;
    for(int i=0;i<m*n;i++){
        res[i].real = z.real;
        res[i].imaginary = z.imaginary;
        z.real++;
        z.imaginary++;
    }
    return res;
}

void print_matrice_complexe_float_t(complexe_float_t* a, int m, int n){
    printf("Matrice :\n");
    for(int i=0;i<m*n;i++){
        if((i%n==0) && (i!=0)){
            printf("\n");
        }
        printf("(%f,%fi)  ",a[i].real,a[i].imaginary);
    }
    printf("\n");
}

complexe_double_t* matrice_complexe_double_t(int m, int n){
    complexe_double_t* res = malloc(m*n*sizeof(complexe_double_t));
    complexe_double_t z;
    z.real = 1.0;
    z.imaginary = 1.0;
    for(int i=0;i<m*n;i++){
        res[i].real = z.real;
        res[i].imaginary = z.imaginary;
        z.real++;
        z.imaginary++;
    }
    return res;
}

void print_matrice_complexe_double_t(complexe_double_t* a, int m, int n){
    printf("Matrice :\n");
    for(int i=0;i<m*n;i++){
        if((i%n==0) && (i!=0)){
            printf("\n");
        }
        printf("(%lf,%lfi)  ",a[i].real,a[i].imaginary);
    }
    printf("\n");
}


int main (int argc, char **argv){

    //MATRICE FLOAT

    float* A = matrice_float(2,2);
    printf("A : ");
    print_matrice_float(A,2,2);
    float* B = matrice_float(2,2);
    printf("B : ");
    print_matrice_float(B,2,2);
    float* C = matrice_float(2,2);
    printf("C : ");
    print_matrice_float(C,2,2);
    float alpha = 2;
    float beta = 2;
    /*float* tmp = MultiplyByScalarConstMat(alpha,A,2,2);
    printf("tmp : ");
    print_matrice_float(tmp,2,2);
    float* tmp2 = ProduitMatriciel(tmp,B,2,2,2);
    printf("produit :");
    print_matrice_float(tmp2,2,2);
    MultiplyByScalar(beta,C,2,2);
    printf("C : ");
    print_matrice_float(C,2,2);
    AdditionMatriciel(tmp2,C,2,2);
    printf("res :");
    print_matrice_float(C,2,2);*/
    mncblas_sgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,2,2,2,alpha,A,-1,B,-1,beta,C,-1);
    printf("RES :");
    print_matrice_float(C,2,2);
    free(A);
    free(B);
    free(C);

    //MATRICE DOUBLE

    double* D = matrice_double(2,2);
    printf("D : ");
    print_matrice_double(D,2,2);
    double* E = matrice_double(2,2);
    printf("E : ");
    print_matrice_double(E,2,2);
    double* F = matrice_double(2,2);
    printf("F : ");
    print_matrice_double(F,2,2);
    double alphad = 2;
    double betad = 2;
    mncblas_dgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,2,2,2,alphad,D,-1,E,-1,betad,F,-1);
    printf("RES :");
    print_matrice_double(F,2,2);
    free(D);
    free(E);
    free(F);

    //MATRICE COMPLEXE FLOAT

    complexe_float_t* G = matrice_complexe_float_t(2,2);
    printf("G : ");
    print_matrice_complexe_float_t(G,2,2);
    complexe_float_t* H = matrice_complexe_float_t(2,2);
    printf("H : ");
    print_matrice_complexe_float_t(H,2,2);
    complexe_float_t* I = matrice_complexe_float_t(2,2);
    printf("I : ");
    print_matrice_complexe_float_t(I,2,2);
    complexe_float_t alpha_complexe;
    alpha_complexe.real = 2.0;
    alpha_complexe.imaginary = 2.0;
    complexe_float_t beta_complexe;
    beta_complexe.real = 2.0;
    beta_complexe.imaginary = 2.0;
    void* alpha_complexe_cast = &alpha_complexe;
    void* beta_complexe_cast = &beta_complexe;

    mncblas_cgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,2,2,2,alpha_complexe_cast,G,-1,H,-1,beta_complexe_cast,I,-1);
    printf("RES :");
    print_matrice_complexe_float_t(I,2,2);
    free(G);
    free(H);
    free(I);

    //MATRICE COMPLEXE DOUBLE

    complexe_double_t* J = matrice_complexe_double_t(2,2);
    printf("J : ");
    print_matrice_complexe_double_t(J,2,2);
    complexe_double_t* K = matrice_complexe_double_t(2,2);
    printf("K : ");
    print_matrice_complexe_double_t(K,2,2);
    complexe_double_t* L = matrice_complexe_double_t(2,2);
    printf("L : ");
    print_matrice_complexe_double_t(L,2,2);
    complexe_double_t alpha_complexe_double;
    alpha_complexe_double.real = 2.0;
    alpha_complexe_double.imaginary = 2.0;
    complexe_double_t beta_complexe_double;
    beta_complexe_double.real = 2.0;
    beta_complexe_double.imaginary = 2.0;
    void* alpha_complexe_double_cast = &alpha_complexe_double;
    void* beta_complexe_double_cast = &beta_complexe_double;

    mncblas_zgemm(MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,2,2,2,alpha_complexe_double_cast,J,-1,K,-1,beta_complexe_double_cast,L,-1);
    printf("RES :");
    print_matrice_complexe_double_t(L,2,2);
    free(J);
    free(K);
    free(L);

}