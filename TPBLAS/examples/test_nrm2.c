#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

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

    float f1[3];
    f1[0] = 4.3;
    f1[1] = 10.2;
    f1[2] = 2.8;
    printf("f1\n");
    print_vectf(f1,3);

    float res;
    unsigned long long int start, end ;

    start =_rdtsc () ;
 
    for (int i = 0 ; i < NB_FOIS; i++){
        res = mnblas_snrm2(3,f1,1);
        f1[1] +=0.1;
    }
    //res = mnblas_snrm2(3,f1,1);

    end = _rdtsc () ;

    printf ("apres boucle res : %f %lld cycles \n", res, end-start) ;

    calcul_flop ("nrm2 float ", NB_FOIS*4, end-start) ;

    double d1[3];
    d1[0] = 4.3;
    d1[1] = 10.2;
    d1[2] = 2.8;

    printf("d1\n");
    print_vectd(d1,3);

    double resd;
    start =_rdtsc () ;

    for (int i = 0 ; i < NB_FOIS; i++){
        resd = mnblas_dnrm2(3,d1,1);
        d1[2] +=0.1;
    }
    //resd = mnblas_dnrm2(3,d1,1);

    end = _rdtsc () ;

    printf ("apres boucle resd : %lf %lld cycles \n", resd, end-start) ;

    calcul_flop ("nrm2 double ", NB_FOIS*4, end-start) ;


    complexe_float_t c1[3];
    c1[0].real = 4.3;
    c1[0].imaginary = 1.2;
    c1[1].real = 10.2;
    c1[1].imaginary = 3.4;
    c1[2].real = 2.8;
    c1[2].imaginary = 14.23;

    printf("c1\n");
    print_vectcf(c1,3);

    float resc;
    start =_rdtsc () ;

    for (int i = 0 ; i < NB_FOIS; i++){
        resc = mnblas_scnrm2(3,c1,1);
        c1[1].real +=0.1;
    }
    //resc = mnblas_scnrm2(3,c1,1);
    end = _rdtsc () ;

    printf ("apres boucle resc : %f %lld cycles \n", resc, end-start) ;

    calcul_flop ("nrm2 complexe float ", NB_FOIS*6, end-start) ;

    complexe_double_t z1[3];
    z1[0].real = 4.3;
    z1[0].imaginary = 1.2;
    z1[1].real = 10.2;
    z1[1].imaginary = 3.4;
    z1[2].real = 2.8;
    z1[2].imaginary = 14.23;

    printf("z1\n");
    print_vectcd(z1,3);

    double rescd;
    start =_rdtsc () ;

    for (int i = 0 ; i < NB_FOIS; i++){
        rescd =  mnblas_dznrm2(3,z1,1);
        z1[0].real += 0.1;
    }
    //rescd = mnblas_dznrm2(3,z1,1);
    end = _rdtsc () ;

    printf ("apres boucle rescd : %lf %lld cycles \n", rescd, end-start) ;

    calcul_flop ("nrm2 complexe double ", NB_FOIS*6, end-start) ;

    exit (0) ;
}
