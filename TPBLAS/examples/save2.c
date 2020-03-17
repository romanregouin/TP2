#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        1000

void AfficherVecteurFloat (int taille, float* vect) {
    for(int i = 0; i != taille; i++) printf("%.2f\n", vect[i]);
    printf("\n");
}

void AfficherVecteurDouble (int taille, double* vect) {
    for(int i = 0; i != taille; i++) printf("%.2f\n", vect[i]);
    printf("\n");
}

void AfficherVecteurFloatC (int taille, complexe_float_t* vect) {
    for(int i = 0; i != taille; i++) printf("%.2f + i*%2.f\n", vect[i].real, vect[i].imaginary);
    printf("\n");
}

void AfficherVecteurDoubleC (int taille, complexe_double_t* vect) {
    for(int i = 0; i != taille; i++) printf("%.2f + i*%2.f\n", vect[i].real, vect[i].imaginary);
    printf("\n");
}

int main (int argc, char **argv) {
    init_flop ();
    unsigned long long int start, end;
    int incX = 1;
    int incY = 1;
    int M = 36;
    int N = 36;



    printf("\t\t\t VECTEUR FLOAT \n");

    float Af[M * N];
    for(int i = 0; i != M * N; i++) Af[i] = i % 17;
    float Xf[M];
    for(int i = 0; i != M; i++) Xf[i] = 1;
    float Yf[M];
    for(int i = 0; i != M; i++) Yf[i] = 1;
    float alphaf = 2;
    float betaf = 0.5;

    start =_rdtsc ();
    mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, alphaf, Af, M, Xf, incX, betaf, Yf, incY);
    end = _rdtsc ();
    //printf("Row-major, NoTrans :\n");
    //AfficherVecteurFloat(M, Yf);
    printf ("sgemv (Row-major, NoTrans) avec %lld cycles \n", end - start);
    calcul_flop ("sgemv (Row-major, NoTrans) :", 2 * M + 7 * M * N + M, end-start);

    for(int i = 0; i != M; i++) Yf[i] = 1;

    start =_rdtsc ();
    mncblas_sgemv(MNCblasRowMajor, MNCblasTrans, M, N, alphaf, Af, M, Xf, incX, betaf, Yf, incY);
    end = _rdtsc ();
    //printf("Row-major, Trans :\n");
    //AfficherVecteurFloat(M, Yf);
    printf ("\nsgemv (Row-major, Trans) avec %lld cycles \n", end - start);
    calcul_flop ("sgemv (Row-major, Trans) :", 2 * M + 7 * M * N + M, end-start);

    for(int i = 0; i != M; i++) Yf[i] = 1;

    start =_rdtsc ();
    mncblas_sgemv(MNCblasRowMajor, MNCblasConjTrans, M, N, alphaf, Af, M, Xf, incX, betaf, Yf, incY);
    end = _rdtsc ();
    //printf("Row-major, ConjTrans :\n");
    //AfficherVecteurFloat(M, Yf);
    printf ("\nsgemv (Row-major, ConjTrans) avec %lld cycles \n", end - start);
    calcul_flop ("sgemv (Row-major, ConjTrans) :", 2 * M + 7 * M * N + M, end-start);



    printf("\n\n\n\t\t\t VECTEUR DOUBLE \n");

    double Ad[M * N];
    for(int i = 0; i != M * N; i++) Ad[i] = i % 17;
    double Xd[M];
    for(int i = 0; i != M; i++) Xd[i] = 1;
    double Yd[M];
    for(int i = 0; i != M; i++) Yd[i] = 1;
    double alphad = 2;
    double betad = 0.5;

    start =_rdtsc ();
    mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, alphad, Ad, M, Xd, incX, betad, Yd, incY);
    end = _rdtsc ();
    //printf("Row-major, NoTrans :\n");
    //AfficherVecteurDouble(M, Yd);
    printf ("dgemv (Row-major, NoTrans) avec %lld cycles \n", end - start);
    calcul_flop ("dgemv (Row-major, NoTrans) :", 2 * M + 7 * M * N + M, end-start);

    for(int i = 0; i != M; i++) Yd[i] = 1;

    start =_rdtsc ();
    mncblas_dgemv(MNCblasRowMajor, MNCblasTrans, M, N, alphad, Ad, M, Xd, incX, betad, Yd, incY);
    end = _rdtsc ();
    //printf("Row-major, Trans :\n");
    //AfficherVecteurDouble(M, Yd);
    printf ("\ndgemv (Row-major, Trans) avec %lld cycles \n", end - start);
    calcul_flop ("dgemv (Row-major, Trans) :", 2 * M + 7 * M * N + M, end-start);

    for(int i = 0; i != M; i++) Yd[i] = 1;

    start =_rdtsc ();
    mncblas_dgemv(MNCblasRowMajor, MNCblasConjTrans, M, N, alphad, Ad, M, Xd, incX, betad, Yd, incY);
    end = _rdtsc ();
    //printf("Row-major, ConjTrans :\n");
    //AfficherVecteurDouble(M, Yd);
    printf ("\ndgemv (Row-major, ConjTrans) avec %lld cycles \n", end - start);
    calcul_flop ("dgemv (Row-major, ConjTrans) :", 2 * M + 7 * M * N + M, end-start);



    printf("\n\n\n\t\t\t VECTEUR COMPLEXE FLOAT \n");

    complexe_float_t Afc[36 * 36];
    for(int i = 0; i != 36 * 36; i++) {
        Afc[i].real = i % 17;
        Afc[i].imaginary = 0;
    }
    complexe_float_t Xfc[36];
    for(int i = 0; i != 36; i++) {
        Xfc[i].real = 1;
        Xfc[i].imaginary = 0;
    }
    complexe_float_t Yfc[36];
    for(int i = 0; i != 36; i++) {
        Yfc[i].real = 1;
        Yfc[i].imaginary = 0;
    }
    complexe_float_t alphafc[1];
    alphafc[0].real = 2;
    alphafc[0].imaginary = 0;
    complexe_float_t betafc[1];
    betafc[0].real = 0.5;
    betafc[0].imaginary = 0;


    start =_rdtsc ();
    mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, alphafc, Afc, M, Xfc, incX, betafc, Yfc, incY);
    end = _rdtsc ();
    //printf("Row-major, NoTrans :\n");
    //AfficherVecteurFloatC(M, Yfc);
    printf ("cgemv (Row-major, NoTrans) avec %lld cycles \n", end - start);
    calcul_flop ("cgemv (Row-major, NoTrans) :", 8 * M + (10 + 7 + 2) * N * M + M * 2, end-start);

    for(int i = 0; i != M; i++) {
        Yfc[i].real = 1;
        Yfc[i].imaginary = 0;
    }

    start =_rdtsc ();
    mncblas_cgemv(MNCblasRowMajor, MNCblasTrans, M, N, alphafc, Afc, M, Xfc, incX, betafc, Yfc, incY);
    end = _rdtsc ();
    //printf("Row-major, Trans :\n");
    //AfficherVecteurFloatC(M, Yfc);
    printf ("\ncgemv (Row-major, Trans) avec %lld cycles \n", end - start);
    calcul_flop ("cgemv (Row-major, Trans) :", 8 * N + (10 + 7 + 2) * M * N + N * 2, end-start);

    for(int i = 0; i != M; i++) {
        Yfc[i].real = 1;
        Yfc[i].imaginary = 0;
    }

    start =_rdtsc ();
    mncblas_cgemv(MNCblasRowMajor, MNCblasConjTrans, M, N, alphafc, Afc, M, Xfc, incX, betafc, Yfc, incY);
    end = _rdtsc ();
    //printf("Row-major, ConjTrans :\n");
    //AfficherVecteurFloatC(M, Yfc);
    printf ("\ncgemv (Row-major, ConjTrans) avec %lld cycles \n", end - start);
    calcul_flop ("cgemv (Row-major, ConjTrans) :", 8 * N + (2 + 3 + 8 + 7 + 2) * M * N + N * 2, end-start);



    printf("\n\n\n\t\t\t VECTEUR COMPLEXE DOUBLE \n");

    complexe_float_t Adc[M * N];
    for(int i = 0; i != M * N; i++) {
        Adc[i].real = i % 17;
        Adc[i].imaginary = 0;
    }
    complexe_float_t Xdc[M];
    for(int i = 0; i != M; i++) {
        Xdc[i].real = 1;
        Xdc[i].imaginary = 0;
    }
    complexe_float_t Ydc[M];
    for(int i = 0; i != M; i++) {
        Ydc[i].real = 1;
        Ydc[i].imaginary = 0;
    }
    complexe_float_t alphadc[1];
    alphadc[0].real = 2;
    alphadc[0].imaginary = 0;
    complexe_float_t betadc[1];
    betadc[0].real = 0.5;
    betadc[0].imaginary = 0;


    start =_rdtsc ();
    mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, alphadc, Adc, M, Xdc, incX, betadc, Ydc, incY);
    end = _rdtsc ();
    //printf("Row-major, NoTrans :\n");
    //AfficherVecteurDoubleC(M, Ydc);
    printf ("zgemv (Row-major, NoTrans) avec %lld cycles \n", end - start);
    calcul_flop ("zgemv (Row-major, NoTrans) :", 8 * M + (10 + 7 + 2) * N * M + M * 2, end-start);

    for(int i = 0; i != M; i++) {
        Ydc[i].real = 1;
        Ydc[i].imaginary = 0;
    }

    start =_rdtsc ();
    mncblas_zgemv(MNCblasRowMajor, MNCblasTrans, M, N, alphadc, Adc, M, Xdc, incX, betadc, Ydc, incY);
    end = _rdtsc ();
    //printf("Row-major, Trans :\n");
    //AfficherVecteurDoubleC(M, Ydc);
    printf ("\nzgemv (Row-major, Trans) avec %lld cycles \n", end - start);
    calcul_flop ("zgemv (Row-major, Trans) :", 8 * N + (10 + 7 + 2) * M * N + N * 2, end-start);

    for(int i = 0; i != M; i++) {
        Ydc[i].real = 1;
        Ydc[i].imaginary = 0;
    }

    start =_rdtsc ();
    mncblas_zgemv(MNCblasRowMajor, MNCblasConjTrans, M, N, alphadc, Adc, M, Xdc, incX, betadc, Ydc, incY);
    end = _rdtsc ();
    //printf("Row-major, ConjTrans :\n");
    //AfficherVecteurDoubleC(M, Ydc);
    printf ("\nzgemv (Row-major, ConjTrans) avec %lld cycles \n", end - start);
    calcul_flop ("zgemv (Row-major, ConjTrans) :", 8 * N + (2 + 3 + 8 + 7 + 2) * M * N + N * 2, end-start);

    return 0;
}