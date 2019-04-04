/*  Copyright (c) 2016-2017 Drew Schmidt
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

#ifndef __COOP_LIB_INVERSE_H__
#define __COOP_LIB_INVERSE_H__


#include <stdlib.h>

#include "lapack.h"
#include "safeomp.h"


static inline int inv_gen_lu(const int n, double *x)
{
  int info = 0;
  
  // Factor x = LU
  int *ipiv = malloc(n * sizeof(*ipiv));
  dgetrf_(&n, &n, x, &n, ipiv, &info);
  if (info != 0) goto cleanup;
  
  // Invert
  int lwork = -1;
  double lwork_dbl;
  dgetri_(&n, x, &n, ipiv, &lwork_dbl, &lwork, &info);
  if (info != 0) goto cleanup;
  
  lwork = (int) lwork_dbl;
  double *work = malloc(lwork * sizeof(*work));
  dgetri_(&n, x, &n, ipiv, work, &lwork, &info);
  
  
  free(work);
  cleanup:
  free(ipiv);
  
  return info;
}



// invert triangular x in place
static inline int inv_tri(const int lda, const int n, double *x)
{
  int info = 0; // initial to shut up static analyzer
  const char uplo = 'u';
  const char diag = 'n';
  
  dtrtri_(&uplo, &diag, &n, x, &lda, &info);
  return info;
}



// invert diagonal x and square it in place
static inline void inv_diagsq(const int len, double *x)
{
  SAFE_FOR_SIMD
  for (int i=0; i<len; i++)
  {
    const double tmp = x[i];
    x[i] = 1.0 / (tmp * tmp);
  }
}



// invert symmetric positive definite (e.g., covariance matrix) via cholesky
static inline int inv_sym_chol(const int n, double *x)
{
  int info;
  const char uplo = 'l';
  
  // factor 
  dpotrf_(&uplo, &n, x, &n, &info);
  if (info != 0)
    return info;
  
  // invert
  dpotri_(&uplo, &n, x, &n, &info);
  
  return info;
}



#endif
