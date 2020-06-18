#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    1000

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

 init_flop () ;
 float res1 ;
 vectorf_init (vec_f1, 1.0) ;
 vectorf_init (vec_f2, 2.0) ;
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res1 = 0.0 ;
     
     start = _rdtsc () ;
        res1 = mncblas_sdot (VECSIZE, vec_f1, 1, vec_f2, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res1, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }
  double res2 ;
  vectord_init (vec_d1, 1.0) ;
  vectord_init (vec_d2, 2.0) ;
  for (i = 0 ; i < NB_FOIS; i++)
   {
     
     res2 = 0.0 ;
     
     start = _rdtsc () ;
        res2 = mncblas_ddot (VECSIZE, vec_d1, 1, vec_d2, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_ddot %d : res = %3.2f nombre de cycles: %Ld \n", i, res2, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }


  complexe_float_t res3;
  complexe_double_t res4;
  vectorcf_init (vec_cf1, 1.0) ;
  vectorcf_init (vec_cf2, 2.0) ;
  for (i = 0 ; i < NB_FOIS; i++)
   {
     
     res3.imaginary = 0.0 ;
     res3.real = 0.0;
     start = _rdtsc () ;
         mncblas_cdotu_sub (VECSIZE, vec_cf1, 1, vec_cf2, 1,&res3) ;
     end = _rdtsc () ;
     
     printf ("mncblas_cdotu_sub %d : res.real = %3.2f res.im = %3.2f nombre de cycles: %Ld \n", i, res3.real,res3.imaginary, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }

   for (i = 0 ; i < NB_FOIS; i++)
   {
     
     res3.imaginary = 0.0 ;
     res3.real = 0.0;
     start = _rdtsc () ;
         mncblas_cdotc_sub (VECSIZE, vec_cf1, 1, vec_cf2, 1,&res3) ;
     end = _rdtsc () ;
     
     printf ("mncblas_cdotc_sub %d : res.real = %3.2f res.im = %3.2f nombre de cycles: %Ld \n", i, res3.real,res3.imaginary, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }

  vectorcd_init (vec_cd1, 1.0) ;
  vectorcd_init (vec_cd2, 2.0) ;
  for (i = 0 ; i < NB_FOIS; i++)
   {
     res4.imaginary = 0.0 ;
     res4.real = 0.0;
     
     
     
     start = _rdtsc () ;
         mncblas_zdotu_sub (VECSIZE, vec_cd1, 1, vec_cd2, 1,&res4) ;
     end = _rdtsc () ;
     
     printf ("mncblas_zdotu_sub %d : res.real = %3.2f res.im = %3.2f nombre de cycles: %Ld \n", i, res4.real,res4.imaginary, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }

   for (i = 0 ; i < NB_FOIS; i++)
   {
     res4.imaginary = 0.0 ;
     res4.real = 0.0;
     
     
     
     start = _rdtsc () ;
         mncblas_zdotc_sub (VECSIZE, vec_cd1, 1, vec_cd2, 1,&res4) ;
     end = _rdtsc () ;
     
     printf ("mncblas_zdotc_sub %d : res.real = %3.2f res.im = %3.2f nombre de cycles: %Ld \n", i, res4.real,res4.imaginary, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }
}
