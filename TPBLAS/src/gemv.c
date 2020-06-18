#include <omp.h>

#include "complexe2.h"
#include "mnblas.h"

void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA,
                   const int M, const int N, const float alpha, const float *A,
                   const int lda, const float *X, const int incX,
                   const float beta, float *Y, const int incY) {
  register unsigned int i;
  register unsigned int j;
  int i_m;
  float tmp;
  if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) ||
      (layout == MNCblasColMajor && TransA != MNCblasNoTrans)) {
    for (i = 0; i < M; i += incY) {
      i_m = i * M;
      tmp = 0;
#pragma omp parallel
#pragma omp for
      for (j = 0; j < N; j += 8) {
        tmp += *(A + i_m + j) * (*(X + j));
        tmp += *(A + i_m + j + 1) * (*(X + j + 1));
        tmp += *(A + i_m + j + 2) * (*(X + j + 2));
        tmp += *(A + i_m + j + 3) * (*(X + j + 3));
        tmp += *(A + i_m + j + 4) * (*(X + j + 4));
        tmp += *(A + i_m + j + 5) * (*(X + j + 5));
        tmp += *(A + i_m + j + 6) * (*(X + j + 6));
        tmp += *(A + i_m + j + 7) * (*(X + j + 7));
      }
      *(Y + i) = *(Y + i) * beta + tmp * alpha;
    }
  } else {
    for (i = 0; i < N; i += incY) {
      tmp = 0;
#pragma omp parallel
#pragma omp for
      for (j = 0; j < M; j += 8) {
        tmp += *(A + i + j * M) * (*(X + j));
        tmp += *(A + i + j * M + 1) * (*(X + j + 1));
        tmp += *(A + i + j * M + 2) * (*(X + j + 2));
        tmp += *(A + i + j * M + 3) * (*(X + j + 3));
        tmp += *(A + i + j * M + 4) * (*(X + j + 4));
        tmp += *(A + i + j * M + 5) * (*(X + j + 5));
        tmp += *(A + i + j * M + 6) * (*(X + j + 6));
        tmp += *(A + i + j * M + 7) * (*(X + j + 7));
      }
      *(Y + i) = *(Y + i) * beta + tmp * alpha;
    }
  }
  return;
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M,
                   const int N, const double alpha, const double *A,
                   const int lda, const double *X, const int incX,
                   const double beta, double *Y, const int incY) {
  register unsigned int i;
  register unsigned int j;
  double tmp;
  if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) ||
      (layout == MNCblasColMajor && TransA != MNCblasNoTrans)) {
    for (i = 0; i < M; i += incY) {
      for (j = 0, tmp = 0; (j < N); j += incX) {
        tmp += *(A + i * lda + j) * (*(X + j));
      }
      *(Y + i) = *(Y + i) * beta + tmp * alpha;
    }
  } else {
    for (i = 0; i < N; i += incY) {
      for (j = 0, tmp = 0; (j < M); j += incX) {
        tmp += *(A + i + j * lda) * (*(X + j));
      }
      *(Y + i) = *(Y + i) * beta + tmp * alpha;
    }
  }
  return;
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M,
                   const int N, const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta, void *Y,
                   const int incY) {
  register unsigned int i;
  register unsigned int j;
  complexe_float_t tmp;
  if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) ||
      (layout == MNCblasColMajor && TransA == MNCblasTrans)) {
    for (i = 0; i < M; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < N); j += incX) {
        tmp = add_complexe_float(
            tmp, mult_complexe_float(((complexe_float_t *)A)[i * lda + j],
                                     ((complexe_float_t *)X)[j]));
      }
      ((complexe_float_t *)Y)[i] = add_complexe_float(
          mult_complexe_float(((complexe_float_t *)Y)[i],
                              ((complexe_float_t *)beta)[0]),
          mult_complexe_float(tmp, ((complexe_float_t *)alpha)[0]));
    }
  } else if (layout == MNCblasRowMajor && TransA == MNCblasConjTrans) {
    for (i = 0; i < N; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < M); j += incX) {
        tmp = add_complexe_float(
            tmp, mult_complexe_float(((complexe_float_t *)A)[i + j * lda],
                                     ((complexe_float_t *)X)[j]));
      }
      ((complexe_float_t *)Y)[i] = add_complexe_float(
          mult_complexe_float(((complexe_float_t *)Y)[i],
                              ((complexe_float_t *)beta)[0]),
          mult_complexe_float(tmp, ((complexe_float_t *)alpha)[0]));
    }
  } else if (layout == MNCblasColMajor && TransA == MNCblasConjTrans) {
    for (i = 0; i < M; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < N); j += incX) {
        tmp = add_complexe_float(
            tmp, mult_complexe_float(((complexe_float_t *)A)[i * lda + j],
                                     ((complexe_float_t *)X)[j]));
      }
      ((complexe_float_t *)Y)[i] = add_complexe_float(
          mult_complexe_float(((complexe_float_t *)Y)[i],
                              ((complexe_float_t *)beta)[0]),
          mult_complexe_float(tmp, ((complexe_float_t *)alpha)[0]));
    }
  } else {
    for (i = 0; i < N; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < M); j += incX) {
        tmp = add_complexe_float(
            tmp, mult_complexe_float(((complexe_float_t *)A)[i + j * lda],
                                     ((complexe_float_t *)X)[j]));
      }
      ((complexe_float_t *)Y)[i] = add_complexe_float(
          mult_complexe_float(((complexe_float_t *)Y)[i],
                              ((complexe_float_t *)beta)[0]),
          mult_complexe_float(tmp, ((complexe_float_t *)alpha)[0]));
    }
  }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M,
                   const int N, const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta, void *Y,
                   const int incY) {
  register unsigned int i;
  register unsigned int j;
  complexe_double_t tmp;
  if ((layout == MNCblasRowMajor && TransA == MNCblasNoTrans) ||
      (layout == MNCblasColMajor && TransA == MNCblasTrans)) {
    for (i = 0; i < M; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < N); j += incX) {
        tmp = add_complexe_double(
            tmp, mult_complexe_double(((complexe_double_t *)A)[i * lda + j],
                                      ((complexe_double_t *)X)[j]));
      }
      ((complexe_double_t *)Y)[i] = add_complexe_double(
          mult_complexe_double(((complexe_double_t *)Y)[i],
                               ((complexe_double_t *)beta)[0]),
          mult_complexe_double(tmp, ((complexe_double_t *)alpha)[0]));
    }
  } else if (layout == MNCblasRowMajor && TransA == MNCblasConjTrans) {
    for (i = 0; i < N; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < M); j += incX) {
        tmp = add_complexe_double(
            tmp, mult_complexe_double(((complexe_double_t *)A)[i + j * lda],
                                      ((complexe_double_t *)X)[j]));
      }
      ((complexe_double_t *)Y)[i] = add_complexe_double(
          mult_complexe_double(((complexe_double_t *)Y)[i],
                               ((complexe_double_t *)beta)[0]),
          mult_complexe_double(tmp, ((complexe_double_t *)alpha)[0]));
    }
  } else if (layout == MNCblasColMajor && TransA == MNCblasConjTrans) {
    for (i = 0; i < M; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < N); j += incX) {
        tmp = add_complexe_double(
            tmp, mult_complexe_double(((complexe_double_t *)A)[i * lda + j],
                                      ((complexe_double_t *)X)[j]));
      }
      ((complexe_double_t *)Y)[i] = add_complexe_double(
          mult_complexe_double(((complexe_double_t *)Y)[i],
                               ((complexe_double_t *)beta)[0]),
          mult_complexe_double(tmp, ((complexe_double_t *)alpha)[0]));
    }
  } else {
    for (i = 0; i < N; i += incY) {
      for (j = 0, tmp.real = 0, tmp.imaginary = 0; (j < M); j += incX) {
        tmp = add_complexe_double(
            tmp, mult_complexe_double(((complexe_double_t *)A)[i + j * lda],
                                      ((complexe_double_t *)X)[j]));
      }
      ((complexe_double_t *)Y)[i] = add_complexe_double(
          mult_complexe_double(((complexe_double_t *)Y)[i],
                               ((complexe_double_t *)beta)[0]),
          mult_complexe_double(tmp, ((complexe_double_t *)alpha)[0]));
    }
  }
}