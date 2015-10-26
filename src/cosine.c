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



// replaces the crossproduct of a matrix with its cosine similarity
static void fill(const unsigned int n, double *restrict crossprod)
{
  int i, j;
  double *diag = malloc(n * sizeof(*diag));
  double diagj;
  
  // store the diagonal to avoid cache thrashing
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
    diag[i] = crossprod[i + n*i];
  
  // Fill lower triangle and diagonal
  #pragma omp parallel for private(i,j) default(shared) if(n>OMP_MIN_SIZE) schedule(dynamic)
  for (j=0; j<n; j++)
  {
    diagj = diag[j];
    
    crossprod[j + n*j] = 1.0;
    
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      crossprod[i + n*j] /= sqrt(diag[i] * diagj);
  }
  
  free(diag);
}



// Copy upper triangle to lower
static inline void symmetrize(const int n, double *restrict x)
{
  int i, j;
  
  SAFE_FOR_SIMD
  for (j=0; j<n; j++)
  {
    for (i=j+1; i<n; i++)
      x[j + n*i] = x[i + n*j];
  }
}



// x is mxn, cos is nxn
void cosim(const int m, const int n, double *restrict x, double *restrict cos)
{
  crossprod(m, n, x, 1.0, cos);
  fill(n, cos);
  symmetrize(n, cos);
}

