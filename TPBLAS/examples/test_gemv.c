#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#include "flop.h"

#define VECSIZE    36//65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef float mfloat [VECSIZE*VECSIZE] ;
typedef double vdouble [VECSIZE] ;
typedef complexe_float_t vcfloat [VECSIZE] ;
typedef complexe_double_t vcdouble [VECSIZE] ;
typedef complexe_double_t mcdouble [VECSIZE*VECSIZE];

vfloat vec_f1, vec_f2 ;
mfloat m_f1;
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

void matricef_init (mcdouble V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE*VECSIZE; i++){
    V[i].real = i%17;
    V[i].imaginary=(i+2)%17;
 }
  return ;
}

void vectorcf_init (vcfloat V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    V [i].imaginary = x ;
    V [i].real = x ;
  }
  return ;
}

void vectord_init (vdouble V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vectorcd_init (vcdouble V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    V [i].imaginary = x ;
    V [i].real = x ;
  }
  return ;
}

void vectorf_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vectorcf_print (vcfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f,%f ", V[i].real,V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}

void vectord_print (vdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vectorcd_print (vcdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("re%f,im%f \n",V[i].real,V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}





int main (int argc, char **argv)
{
 unsigned long long start, end ;
 
 int i ;
 float alpha1=2.0;
 float alpha2=0.5;
 complexe_double_t alpha3;
 alpha3.imaginary=0;
 alpha3.real=2;
 complexe_double_t alpha4;
 alpha4.imaginary=0;
 alpha4.real=0.5;
 init_flop () ;
 vectorcd_init (vec_cd1, 1.0) ;
 vectorcd_init (vec_cd2, 1.0) ;
 matricef_init (m_cd1, 1.0) ;
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vectorcd_init (vec_cd2, 1.0) ;

     printf ("Y avant      %d : Y[0] = %3.2f,Y[655355] = %3.2f\n", i,vec_f2[0],vec_f2[VECSIZE-1]) ;

     start = _rdtsc () ;
        mncblas_zgemv (MNCblasColMajor,MNCblasConjTrans,VECSIZE,VECSIZE,&alpha3,m_cd1,VECSIZE,vec_cd1, 1,&alpha4, vec_cd2, 1) ;
     end = _rdtsc () ;
     vectorcd_print(vec_cd2);
     printf ("mnblas_sgemv %d : Y[0] = %3.2f,Y[655355] = %3.2f nombre de cycles: %Ld \n", i,vec_cd2[0].real,vec_cd2[0].imaginary, end-start) ;
     calcul_flop ("sdot ", 28*VECSIZE*VECSIZE+32*VECSIZE, end-start) ;
   }
}