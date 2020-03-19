#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#define    NB_FOIS        4194304

void print_vectf(float* v, float* u, int N){
    for(int i=0;i<N;i++){
        printf("%d : %f   ||  ",i,v[i]);
        printf("%d : %f\n",i,u[i]);
    }
}

void print_vectd(double* v, double* u, int N){
    for(int i=0;i<N;i++){
        printf("%d : %lf   ||  ",i,v[i]);
        printf("%d : %lf\n",i,u[i]);
    }
}

void print_vectcf(complexe_float_t* v, complexe_float_t* u, int N){
    for(int i=0;i<N;i++){
        printf("%d : Real : %f Img : %f   ||  ",i,v[i].real,v[i].imaginary);
        printf("%d : Real : %f Img : %f\n",i,u[i].real,u[i].imaginary);
    }
}

void print_vectcd(complexe_double_t* v, complexe_double_t* u, int N){
    for(int i=0;i<N;i++){
        printf("%d : Real : %f Img : %lf   ||  ",i,v[i].real,v[i].imaginary);
        printf("%d : Real : %f Img : %lf\n",i,u[i].real,u[i].imaginary);
    }
}

int main (int argc, char **argv){
    float incrementf = 0.9;
    float f1[10];
    for(int i=0;i<10;i++){
        f1[i] = 0.9 + incrementf;
        incrementf += 0.1;
    }
    float f2[10];
    for(int i=0;i<10;i++){
        f2[i] = 0.1 + incrementf;
        incrementf += 0.1;
    }

    printf("f1 || f2\n");
    print_vectf(f1,f2,10);
    printf("swaped\n");

    mncblas_sswap(10,f1,1,f2,1);

    printf("f1 || f2\n");
    print_vectf(f1,f2,10);

    double incrementd = 0.9;
    double d1[10];
    for(int i=0;i<10;i++){
        d1[i] = 0.9 + incrementd;
        incrementd += 0.1;
    }
    double d2[10];
    for(int i=0;i<10;i++){
        d2[i] = 0.1 + incrementd;
        incrementd += 0.1;
    }

    printf("d1 || d2\n");
    print_vectd(d1,d2,10);
    printf("swaped\n");

    mncblas_dswap(10,d1,1,d2,1);

    printf("d1 || d2\n");
    print_vectd(d1,d2,10);

    float incrementcf = 0.9;
    complexe_float_t c1[10];
    for(int i=0;i<10;i++){
        c1[i].real = 0.9 + incrementcf;
        c1[i].imaginary = 0.9 + incrementcf;
        incrementcf += 0.1;
    }
    complexe_float_t c2[10];
    for(int i=0;i<10;i++){
        c2[i].real = 0.1 + incrementcf;
        c2[i].imaginary = 0.1 + incrementcf;
        incrementcf += 0.1;
    }

    printf("c1 || c2\n");
    print_vectcf(c1,c2,10);
    printf("swaped\n");

    mncblas_cswap(10,c1,1,c2,1);

    printf("c1 || c2\n");
    print_vectcf(c1,c2,10);

    double incrementcd = 0.9;
    complexe_double_t z1[10];
    for(int i=0;i<10;i++){
        z1[i].real = 0.9 + incrementcd;
        z1[i].imaginary = 0.9 + incrementcd;
        incrementcd += 0.1;
    }
    complexe_double_t z2[10];
    for(int i=0;i<10;i++){
        z2[i].real = 0.1 + incrementcd;
        z2[i].imaginary = 0.1 + incrementcd;
        incrementcd += 0.1;
    }

    printf("z1 || z2\n");
    print_vectcd(z1,z2,10);
    printf("swaped\n");

    mncblas_zswap(10,z1,1,z2,1);

    printf("z1 || z2\n");
    print_vectcd(z1,z2,10);

    exit (0) ;
}
