#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#define    NB_FOIS        4194304

#include "flop.h"

void print_vectf(float* v, int N){
    for(int i=0;i<N;i++){
        printf("%d : %f\n",i,v[i]);
    }
}

void print_vectd(double* v, int N){
    for(int i=0;i<N;i++){
        printf("%d : %lf\n",i,v[i]);
    }
}

void print_vectcf(complexe_float_t* v, int N){
    for(int i=0;i<N;i++){
        printf("%d : Real : %f Img : %f\n",i,v[i].real,v[i].imaginary);
    }
}

void print_vectcd(complexe_double_t* v, int N){
    for(int i=0;i<N;i++){
        printf("%d : Real : %f Img : %lf\n",i,v[i].real,v[i].imaginary);
    }
}

int main (int argc, char **argv){

    float incrementf = 0.9;
    float f1[10];
    for(int i=0;i<10;i++){
        f1[i] = 0.9 + incrementf;
        incrementf += 0.1;
    }
    
    printf("f1\n");
    print_vectf(f1,10);

    CBLAS_INDEX indice;
    unsigned long long int start, end ;

    start =_rdtsc () ;
 
    for (int i = 0 ; i < NB_FOIS; i++){
        indice = mnblas_isamax(10,f1,1);
        f1[9] +=0.1;
    }

    end = _rdtsc () ;

    printf ("apres boucle max : %ld %lld cycles \n", indice, end-start) ;

    calcul_flop ("iamax float ", NB_FOIS*2, end-start) ;
    
    double incrementd = 0.9;
    double d1[10];
    for(int i=0;i<10;i++){
        d1[i] = 0.9 + incrementd;
        incrementd += 0.1;
    }

    printf("d1\n");
    print_vectd(d1,10);

    CBLAS_INDEX indiced;
    start =_rdtsc () ;

        for (int i = 0 ; i < NB_FOIS; i++){
        indiced = mnblas_idamax(10,d1,1);
        d1[5] +=0.1;
    }

    end = _rdtsc () ;

    printf ("apres boucle indiced : %ld %lld cycles \n", indiced, end-start) ;

    calcul_flop ("iamax double ", NB_FOIS*2, end-start) ;


    float incrementcf = 0.9;
    complexe_float_t c1[10];
    for(int i=0;i<10;i++){
        c1[i].real = 0.9 + incrementcf;
        c1[i].imaginary = 0.9 + incrementcf;
        incrementcf += 0.1;
    }

    printf("c1\n");
    print_vectcf(c1,10);

    CBLAS_INDEX indicecf;
    start =_rdtsc () ;

        for (int i = 0 ; i < NB_FOIS; i++){
        indicecf = mnblas_icamax(10,c1,1);
        c1[1].real +=0.1;
    }

    end = _rdtsc () ;

    printf ("apres boucle indicecf : %ld %lld cycles \n", indicecf, end-start) ;

    calcul_flop ("iamax complexe float ", NB_FOIS*4, end-start) ;

    double incrementcd = 0.9;
    complexe_double_t z1[10];
    for(int i=0;i<10;i++){
        z1[i].real = 0.9 + incrementcd;
        z1[i].imaginary = 0.9 + incrementcd;
        incrementcd += 0.1;
    }

    printf("z1\n");
    print_vectcd(z1,10);

    CBLAS_INDEX indicecd;
    start =_rdtsc () ;

        for (int i = 0 ; i < NB_FOIS; i++){
        indicecd =  mnblas_izamax(10,z1,1);
        z1[4].real += 0.1;
    }

    end = _rdtsc () ;

    printf ("apres boucle indicecd : %ld %lld cycles \n", indicecd, end-start) ;

    calcul_flop ("iamax complexe double ", NB_FOIS*4, end-start) ;

    exit (0) ;
}
