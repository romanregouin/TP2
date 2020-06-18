#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#include "flop.h"

#define VECSIZE    4//65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef float mfloat [VECSIZE*VECSIZE] ;
typedef double vdouble [VECSIZE] ;
typedef complexe_float_t vcfloat [VECSIZE] ;
typedef complexe_double_t vcdouble [VECSIZE] ;
typedef complexe_float_t mcfloat [VECSIZE*VECSIZE];
typedef complexe_double_t mcdouble [VECSIZE*VECSIZE];

vfloat vec_f1, vec_f2 ;
mfloat m_f1;
mcfloat m_cf1;
mcdouble m_cd1;
vdouble vec_d1, vec_d2 ;
vcfloat vec_cf1, vec_cf2 ;
vcdouble vec_cd1, vec_cd2 ;

void vectorf_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void matricef_init (mcfloat V, float x){
  register unsigned int i ;

  for (i = 0; i < VECSIZE*VECSIZE; i++){
    V[i].real = i%17;
    V[i].imaginary=(i)%17;
 }
  return ;
}


void matriced_init (mcdouble V, double x){
  register unsigned int i ;

  for (i = 0; i < VECSIZE*VECSIZE; i++){
    V[i].real = i%17;
    V[i].imaginary=(i)%17;
 }
  return ;
}


void vectorcf_init (vcfloat V, double x){
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    V [i].imaginary = x ;
    V [i].real = x ;
  }
  return ;
}

void vectord_init (vdouble V, float x){
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vectorcd_init (vcdouble V, double x){
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    V [i].imaginary = x ;
    V [i].real = x ;
  }
  return ;
}

void vectorf_print (vfloat V){
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vectorcf_print (vcfloat V){
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("re%f,im%f \n", V[i].real,V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}

void vectord_print (vdouble V){
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vectorcd_print (vcdouble V){
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("re%f,im%f \n",V[i].real,V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}





int main (int argc, char **argv){
 unsigned long long start, end ;
 
 int i ;
 //float alpha1=2.0;
 //float alpha2=0.5;
 complexe_double_t alpha3;
 alpha3.imaginary=0;
 alpha3.real=2;
 complexe_double_t alpha4;
 alpha4.imaginary=0;
 alpha4.real=0.5;
 init_flop () ;
 vectorcd_init (vec_cd1, 1.0) ;
 vectorcd_init (vec_cd2, 1.0) ;
 matriced_init (m_cd1, 1.0) ;
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vectorcd_init (vec_cd2, 1.0) ;

     

     start = _rdtsc () ;
        mncblas_zgemv (MNCblasRowMajor,MNCblasNoTrans,VECSIZE,VECSIZE,&alpha3,m_cd1,VECSIZE,vec_cd1, 1,&alpha4, vec_cd2, 1) ;
     end = _rdtsc () ;
     vectorcd_print(vec_cd2);
     printf ("mnblas_sgemv %d : Y[0] = %3.2f,Y[655355] = %3.2f nombre de cycles: %Ld \n", i,vec_cd2[0].real,vec_cd2[0].imaginary, end-start) ;
     calcul_flop ("sdot ", 28*VECSIZE*VECSIZE+32*VECSIZE, end-start) ;
   }


   complexe_float_t alpha3_cf;
   alpha3_cf.imaginary=0;
   alpha3_cf.real=2;
   complexe_float_t alpha4_cf;
   alpha4_cf.imaginary=0;
   alpha4_cf.real=0.5;
   init_flop () ;
   vectorcf_init (vec_cf1, 1.0) ;
   vectorcf_init (vec_cf2, 1.0) ;
   matricef_init (m_cf1, 1.0) ;
  for (i = 0 ; i < NB_FOIS; i++)
   {
     vectorcf_init (vec_cf2, 1.0) ;

     
    
     start = _rdtsc () ;
        mncblas_cgemv (MNCblasRowMajor,MNCblasNoTrans,VECSIZE,VECSIZE,&alpha3_cf,m_cf1,VECSIZE,vec_cf1, 1,&alpha4_cf, vec_cf2, 1) ;
     end = _rdtsc () ;
     vectorcf_print(vec_cf2);
     printf ("mnblas_sgemv %d : Y[0] = %3.2f,Y[655355] = %3.2f nombre de cycles: %Ld \n", i,vec_cf2[0].real,vec_cf2[0].imaginary, end-start) ;
     calcul_flop ("sdot ", 28*VECSIZE*VECSIZE+32*VECSIZE, end-start) ;
   }
}