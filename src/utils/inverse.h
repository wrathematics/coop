/*  Copyright (c) 2016, Schmidt
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


#include "lapack.h"
#include "safeomp.h"


static inline int inv_ge_lu(const int n, double *x)
{
  int info = 0;
  int *ipiv;
  int lwork;
  double tmp;
  double *work;
  
  
  // Factor x = LU
  ipiv = malloc(n * sizeof(*ipiv));
  dgetrf_(&n, &n, x, &n, ipiv, &info);
  if (info != 0) goto cleanup;
  
  // Invert
  lwork = -1;
  dgetri_(&n, x, &n, ipiv, &tmp, &lwork, &info);
  if (info != 0) goto cleanup;
  
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  dgetri_(&n, x, &n, ipiv, work, &lwork, &info);
  
  
  free(work);
  cleanup:
  free(ipiv);
  
  return info;
}



// invert triangular x in place
static inline int inv_tr(const int lda, const int n, double *x)
{
  int info;
  const char uplo = 'u';
  const char diag = 'n';
  
  dtrtri_(&uplo, &diag, &n, x, &lda, &info);
  return info;
}



// invert diagonal x in place
static inline void inv_diag(const int len, double *x)
{
  SAFE_FOR_SIMD
  for (int i=0; i<len; i++)
  {
    const register double tmp = x[i];
    x[i] = 1.0 / (tmp * tmp);
  }
}


#endif
