#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#define    NB_FOIS  10000
#define    VECT_SIZE    5000

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
    float f1[VECT_SIZE];
    for(int i=0;i<VECT_SIZE;i++){
        f1[i] = 0.9 + incrementf;
        incrementf += 0.1;
    }
    
    //printf("f1\n");
    //print_vectf(f1,VECT_SIZE);
    printf("summed\n");

    float res;
    unsigned long long int start, end ;

    start =_rdtsc () ;
 
    for (int i = 0 ; i < NB_FOIS; i++){
        res = mnblas_sasum(VECT_SIZE,f1,1);
        f1[0] +=0.1;
        //printf("res : %f\n",res);
    }

    end = _rdtsc () ;

    printf ("apres boucle res : %f %lld cycles \n", res, end-start) ;

    calcul_flop ("asum float ", NB_FOIS*VECT_SIZE, end-start) ;

    

    double incrementd = 0.9;
    double d1[VECT_SIZE];
    for(int i=0;i<VECT_SIZE;i++){
        d1[i] = 0.9 + incrementd;
        incrementd += 0.1;
    }

    //printf("d1\n");
    //print_vectd(d1,VECT_SIZE);
    printf("summed\n");

    double resd;
    start =_rdtsc () ;

        for (int i = 0 ; i < NB_FOIS; i++){
        resd = mnblas_dasum(VECT_SIZE,d1,1);
        d1[0] +=0.1;
        //printf("resd : %lf\n",resd);
    }

    end = _rdtsc () ;

    printf ("apres boucle resd : %f %lld cycles \n", resd, end-start) ;

    calcul_flop ("asum double ", NB_FOIS*VECT_SIZE, end-start) ;

    float incrementcf = 0.9;
    complexe_float_t c1[VECT_SIZE];
    for(int i=0;i<VECT_SIZE;i++){
        c1[i].real = 0.9 + incrementcf;
        c1[i].imaginary = 0.9 + incrementcf;
        incrementcf += 0.1;
    }

    //printf("c1\n");
    //print_vectcf(c1,VECT_SIZE);
    printf("summed\n");

    double resc;
    start =_rdtsc () ;

        for (int i = 0 ; i < NB_FOIS; i++){
        resc = mnblas_scasum(VECT_SIZE,c1,1);
        c1[0].real +=0.1;
        //printf("resc : %f\n",resc);
    }

    end = _rdtsc () ;

    printf ("apres boucle resd : %f %lld cycles \n", resc, end-start) ;

    calcul_flop ("asum complexe float ", NB_FOIS*VECT_SIZE, end-start) ;


    double incrementcd = 0.9;
    complexe_double_t z1[VECT_SIZE];
    for(int i=0;i<VECT_SIZE;i++){
        z1[i].real = 0.9 + incrementcd;
        z1[i].imaginary = 0.9 + incrementcd;
        incrementcd += 0.1;
    }

    //printf("z1\n");
    //print_vectcd(z1,VECT_SIZE);
    printf("summed\n");

    double rescd;
    start =_rdtsc () ;

        for (int i = 0 ; i < NB_FOIS; i++){
        rescd =  mnblas_dzasum(VECT_SIZE,z1,1);
        z1[0].real += 0.1;
        //printf("rescd : %lf\n",rescd);
    }

    end = _rdtsc () ;

    printf ("apres boucle resd : %f %lld cycles \n", rescd, end-start) ;

    calcul_flop ("asum complexe double ", NB_FOIS*VECT_SIZE, end-start) ;
    
    exit (0) ;
}
