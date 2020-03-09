#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef double vdouble [VECSIZE] ;
typedef complexe_float_t vcfloat [VECSIZE] ;
typedef complexe_double_t vcdouble [VECSIZE] ;

vfloat vec_f1, vec_f2 ;
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
    printf ("%f,%f ",V[i].real,V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}


int main (int argc, char **argv)
{
 unsigned long long start, end ;
 
 int i ;
 float alpha1=2.0;
 double alpha2=2.0;
 complexe_float_t alpha3;
 alpha3.imaginary=2.0;
 alpha3.real=0.5;
 complexe_double_t alpha4;
 alpha4.imaginary=2.0;
 alpha4.real=0.5;
 init_flop () ;
 vectorf_init (vec_f1, 1.0) ;
 vectorf_init (vec_f2, 2.0) ;
 for (i = 0 ; i < NB_FOIS; i++)
   {

     printf ("Y avant      %d : Y[0] = %3.2f,Y[655355] = %3.2f\n", i,vec_f2[0],vec_f2[65535]) ;

     start = _rdtsc () ;
        mnblas_saxpy (VECSIZE, alpha1,vec_f1, 1, vec_f2, 1) ;
     end = _rdtsc () ;
     
     printf ("mnblas_saxpy %d : Y[0] = %3.2f,Y[655355] = %3.2f nombre de cycles: %Ld \n", i,vec_f2[0],vec_f2[65535], end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }

  vectord_init (vec_d1, 1.0) ;
  vectord_init (vec_d2, 2.0) ;
  for (i = 0 ; i < NB_FOIS; i++)
   {
     
     printf ("Y avant      %d : Y[0] = %3.2f,Y[655355] = %3.2f\n", i,vec_d2[0],vec_d2[65535]) ;
     
     start = _rdtsc () ;
        mnblas_daxpy (VECSIZE,alpha2, vec_d1, 1, vec_d2, 1) ;
     end = _rdtsc () ;
     
     printf ("mnblas_daxpy %d : Y[0] = %3.2f,Y[655355] = %3.2f nombre de cycles: %Ld \n", i,vec_d2[0],vec_d2[65535], end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }


  vectorcf_init (vec_cf1, 1.0) ;
  vectorcf_init (vec_cf2, 2.0) ;
  for (i = 0 ; i < NB_FOIS; i++)
   {
     
     printf ("Y avant      %d : Y[0].r = %3.2f, Y[0].im = %3.2f  Y[655355].r= %3.2f Y[655355].im= %3.2f\n", i,vec_cf2[0].real,vec_cf2[0].imaginary,vec_cf2[65535].real,vec_cf2[65535].imaginary) ;
     
     start = _rdtsc () ;
        mnblas_caxpy (VECSIZE,&alpha3, vec_cf1, 1, vec_cf2, 1) ;
     end = _rdtsc () ;
     
     printf ("mnblas_caxpy %d : Y[0].r = %3.2f, Y[0].im = %3.2f  Y[655355].r= %3.2f Y[655355].im= %3.2f nombre de cycles: %Ld\n", i,vec_cf2[0].real,vec_cf2[0].imaginary,vec_cf2[65535].real,vec_cf2[65535].imaginary, end-start);
     calcul_flop ("sdot ", 8 * VECSIZE, end-start) ;
   }

  vectorcd_init (vec_cd1, 1.0) ;
  vectorcd_init (vec_cd2, 2.0) ;
  for (i = 0 ; i < NB_FOIS; i++)
   {
     printf ("Y avant      %d : Y[0].r = %3.2f, Y[0].im = %3.2f  Y[655355].r= %3.2f Y[655355].im= %3.2f\n", i,vec_cd2[0].real,vec_cd2[0].imaginary,vec_cd2[65535].real,vec_cd2[65535].imaginary) ;
     
     start = _rdtsc () ;
        mnblas_zaxpy (VECSIZE,&alpha4, vec_cd1, 1, vec_cd2, 1) ;
     end = _rdtsc () ;
     
     printf ("mnblas_zaxpy %d : Y[0].r = %3.2f, Y[0].im = %3.2f  Y[655355].r= %3.2f Y[655355].im= %3.2f nombre de cycles: %Ld\n", i,vec_cd2[0].real,vec_cd2[0].imaginary,vec_cd2[65535].real,vec_cd2[65535].imaginary, end-start);
     calcul_flop ("sdot ", 8 * VECSIZE, end-start) ;
   }

  
}
