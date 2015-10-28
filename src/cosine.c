/*  Copyright (c) 2015, Schmidt
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
    1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <stdlib.h>
#include <math.h>
#include "cosine.h"
#include "omp.h"


void dgemm_(const char *transa, const char *transb, const int *m, const int *n, 
            const int *k, const double *restrict alpha, const double *restrict a, 
            const int *lda, const double *restrict b, const int *ldb, 
            const double *beta, double *restrict c, const int *ldc);


static inline double ddot(const int n, const double *restrict x, const double *restrict y)
{
  int one = 1;
  double dot;
  
  dgemm_(&(char){'t'}, &(char){'n'}, &one, &one, &n, 
    &(double){1.0}, x, &n, y, &n, &(double){0.0}, &dot, &one);
  
  return dot;
}


void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k, 
            const double *restrict alpha, const double *restrict a, const int *lda, 
            const double *restrict beta, double *restrict c, const int *ldc);



// upper triangle of t(x) %*% x
static inline int crossprod(const int m, const int n, const double *restrict x, const double alpha, double *restrict c)
{
  int info = 0;
  
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &alpha, x, &m, &(double){0.0}, c, &n);
  
  return info;
}



// replaces upper triangle of the crossproduct of a matrix with its cosine similarity
static inline void fill(const unsigned int n, double *restrict crossprod)
{
  int i, j;
  double *diag = malloc(n * sizeof(*diag));
  double diagj;
  
  #pragma omp parallel for private(i,j,diagj) default(shared) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    diagj = crossprod[j + n*j];
    
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      crossprod[i + n*j] /= sqrt(crossprod[i + n*i] * diagj);
  }
  
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
    crossprod[i + n*i] = 1.0;
  
  free(diag);
}



// Copy upper triangle to lower
static inline void symmetrize(const int n, double *restrict x)
{
  int i, j;
  
  #pragma omp parallel for private(i) default(shared) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      x[j + n*i] = x[i + n*j];
  }
}



/**
 * @brief 
 * Compute the cosine similarity matrix of a matrix.  This is
 * all pair-wise vector cosine similarities of the columns.
 * 
 * @details
 * The implementation is dominated by a symmetric rank-k update
 * via the BLAS function dsyrk().
 * 
 * @param m,n
 * The number of rows/columns of the input matrix x.
 * @param x
 * The input mxn matrix.
 * @param cos
 * The output nxn matrix.
*/
void cosine_mat(const int m, const int n, const double *restrict x, double *restrict cos)
{
  crossprod(m, n, x, 1.0, cos);
  fill(n, cos);
  symmetrize(n, cos);
}



double cosine_vecvec(const int n, const double *restrict x, const double *restrict y)
{
  double normx, normy;
  
  double cp = ddot(n, x, y);
  
  crossprod(n, 1, x, 1.0, &normx);
  crossprod(n, 1, y, 1.0, &normy);
  
  return cp / sqrt(normx * normy);
}



