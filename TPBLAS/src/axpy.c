#include "mnblas.h"
#include <stdio.h>

void mnblas_saxpy(const int N, const float alpha, const float *X,
                  const int incX, float *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;
    for (; ((i < N) && (j < N)); i += incX, j += incY)
    {
        Y[j] = alpha * X[i] + Y[j];
    }
}

void mnblas_daxpy(const int N, const double alpha, const double *X,
                  const int incX, double *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;
    for (; ((i < N) && (j < N)); i += incX, j += incY)
    {
        Y[j] = alpha * X[i] + Y[j];
    }
}

void mnblas_caxpy(const int N, const void *alpha, const void *X,
                  const int incX, void *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;
    for (; ((i < N) && (j < N)); i += incX, j += incY)
    {
        ((float *)Y)[2*j+1] = ((float *)alpha)[0]*((float *)X)[2*j+1]+((float *)alpha)[1]*((float *)X)[2*j] + ((float *)Y)[2*j+1];
        ((float *)Y)[2*j] = ((float *)alpha)[0]*((float *)X)[2*j]-((float *)alpha)[1]*((float *)X)[2*j+1] + ((float *)Y)[2*j];
      
    }
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X,
                  const int incX, void *Y, const int incY)
{
    register unsigned int i = 0;
    register unsigned int j = 0;
    for (; ((i < N) && (j < N)); i += incX, j += incY)
    {
        ((double *)Y)[2*j+1] = ((double*)alpha)[0]*((double*)X)[2*j+1]+((double*)alpha)[1]*((double*)X)[2*j] + ((double*)Y)[2*j+1];
        ((double *)Y)[2*j] = ((double*)alpha)[0]*((double*)X)[2*j]-((double*)alpha)[1]*((double*)X)[2*j+1] + ((double*)Y)[2*j];
    }
}
