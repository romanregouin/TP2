#include <omp.h>
#include <stdio.h>
#include <x86intrin.h>

#include "complexe2.h"
#include "flop.h"
#include "mnblas.h"

#define VECSIZE 256  // 512

#define NB_FOIS 1

typedef float vfloat[VECSIZE];
typedef float mfloat[VECSIZE * VECSIZE];
typedef double mdouble[VECSIZE * VECSIZE];
typedef double vdouble[VECSIZE];
typedef complexe_float_t vcfloat[VECSIZE];
typedef complexe_double_t vcdouble[VECSIZE];
typedef complexe_float_t mcfloat[VECSIZE * VECSIZE];
typedef complexe_double_t mcdouble[VECSIZE * VECSIZE];

vfloat vec_f1, vec_f2;
mfloat m_f1;
mdouble m_d1;
mcfloat m_cf1;
mcdouble m_cd1;
vdouble vec_d1, vec_d2;
vcfloat vec_cf1, vec_cf2;
vcdouble vec_cd1, vec_cd2;

void vectorf_init(vfloat V, float x) {
  register unsigned int i;

  if (x != -1)
    for (i = 0; i < VECSIZE; i++) V[i] = x - i;
  else
    for (i = 0; i < VECSIZE; i++) V[i] = i;

  return;
}

void matricef_init(mfloat V, float x) {
  register unsigned int i;

  for (i = 0; i < VECSIZE * VECSIZE; i++) {
    V[i] = i;
  }
  return;
}

void matriced_init(mdouble V, float x) {
  register unsigned int i;

  for (i = 0; i < VECSIZE * VECSIZE; i++) {
    V[i] = i;
  }
  return;
}

void matricecf_init(mcfloat V, float x) {
  register signed int i;

  for (i = 0; i < VECSIZE * VECSIZE; i++) {
    V[i].real = i;
    V[i].imaginary = 10 - i;
  }
  return;
}

void matricecd_init(mcdouble V, double x) {
  register signed int i;

  for (i = 0; i < VECSIZE * VECSIZE; i++) {
    V[i].real = i;
    V[i].imaginary = 10 - i;
  }
  return;
}

void vectorcf_init(vcfloat V, double x) {
  register unsigned int i;
  if (x != -1) {
    for (i = 0; i < VECSIZE; i++) {
      V[i].imaginary = i;
      V[i].real = x - i;
    }
  } else {
    for (i = 0; i < VECSIZE; i++) {
      V[i].imaginary = VECSIZE + 1 - i;
      V[i].real = i;
    }
  }
  return;
}

void vectord_init(vdouble V, float x) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++) V[i] = x;

  return;
}

void vectorcd_init(vcdouble V, double x) {
  register unsigned int i;

  if (x != -1) {
    for (i = 0; i < VECSIZE; i++) {
      V[i].imaginary = i;
      V[i].real = x - i;
    }
  } else {
    for (i = 0; i < VECSIZE; i++) {
      V[i].imaginary = VECSIZE + 1 - i;
      V[i].real = i;
    }
  }
  return;
}

void matricecd_print(mcdouble M) {
  register unsigned int i;
  for (i = 0; i < VECSIZE * VECSIZE; i++)
    printf("re %f   im %f   \n ", M[i].real, M[i].imaginary);
  printf("\n");
}

void vectorf_print(vfloat V) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++) printf("%f ", V[i]);
  printf("\n");

  return;
}

void vectorcf_print(vcfloat V) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("re%f,im%f \n", V[i].real, V[i].imaginary);
  printf("\n");

  return;
}

void vectord_print(vdouble V) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++) printf("%f ", V[i]);
  printf("\n");

  return;
}

void vectorcd_print(vcdouble V) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("re%f,im%f \n", V[i].real, V[i].imaginary);
  printf("\n");

  return;
}

int main(int argc, char **argv) {
  unsigned long long start, end;

  int i;
  init_flop();
  printf("%d\n\n\n", omp_get_max_threads());
  vectord_init(vec_d1, 1.0);
  vectord_init(vec_d2, 1.0);
  matriced_init(m_d1, 1.0);
  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++) {
    mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, 2, m_d1,
                  VECSIZE, vec_d1, 1, 0.5, vec_d2, 1);
  }
  end = _rdtsc();
  // vectord_print(vec_d2);
  calcul_flop("dgemv: ", (2 * VECSIZE * VECSIZE + 3 * VECSIZE) * NB_FOIS,
              end - start);
  printf(" %lli\n", end - start);

  vectorf_init(vec_f1, 5);
  vectorf_init(vec_f2, -1);
  matricef_init(m_f1, 1.0);
  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++) {
    mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, 2, m_f1,
                  VECSIZE, vec_f1, 1, 0.5, vec_f2, 1);
  }
  end = _rdtsc();
  // vectorf_print(vec_f2);
  calcul_flop("sgemv (avec deroulement de boucle de 8): ",
              (2 * VECSIZE * VECSIZE + 3 * VECSIZE) * NB_FOIS, end - start);
  printf(" %lli\n", end - start);

  complexe_float_t alpha3_cf;
  alpha3_cf.imaginary = 0;
  alpha3_cf.real = 2;
  complexe_float_t alpha4_cf;
  alpha4_cf.imaginary = 0;
  alpha4_cf.real = 0.5;
  vectorcf_init(vec_cf1, 6);
  vectorcf_init(vec_cf2, -1);
  matricecf_init(m_cf1, 1.0);
  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++) {
    mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, &alpha3_cf,
                  m_cf1, VECSIZE, vec_cf1, 1, &alpha4_cf, vec_cf2, 1);
  }
  end = _rdtsc();
  // vectorcf_print(vec_cf2);
  calcul_flop("cgemv: ", (8 * VECSIZE * VECSIZE + 14 * VECSIZE) * NB_FOIS,
              end - start);
  printf(" %lli\n", end - start);

  complexe_double_t alpha3;
  alpha3.imaginary = 0;
  alpha3.real = 2;
  complexe_double_t alpha4;
  alpha4.imaginary = 0;
  alpha4.real = 0.5;

  vectorcd_init(vec_cd1, 6);
  vectorcd_init(vec_cd2, -1);
  matricecd_init(m_cd1, 1.0);
  // vectorcd_print(vec_cd1);
  // vectorcd_print(vec_cd2);

  start = _rdtsc();
  for (i = 0; i < NB_FOIS; i++) {
    mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, &alpha3,
                  m_cd1, VECSIZE, vec_cd1, 1, &alpha4, vec_cd2, 1);
  }
  end = _rdtsc();
  // vectorcd_print(vec_cd2);
  calcul_flop("zgemv: ", (8 * VECSIZE * VECSIZE + 14 * VECSIZE) * NB_FOIS,
              end - start);
  printf(" %lli\n", end - start);
}