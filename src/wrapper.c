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

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "cosine.h"


SEXP R_cosine_mat(SEXP x)
{
  const unsigned int m = nrows(x);
  const unsigned int n = ncols(x);
  SEXP ret;
  PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  cosine_mat(m, n, REAL(x), REAL(ret));
  
  UNPROTECT(1);
  return ret;
}



SEXP R_cosine_vecvec(SEXP x, SEXP y)
{
  const unsigned int n = LENGTH(x);
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 1));
  
  REAL(ret)[0] = cosine_vecvec(n, REAL(x), REAL(y));
  
  UNPROTECT(1);
  return ret;
}



#define INEDEX_FROM_1 1

SEXP R_cosine_sparse_coo(SEXP n_, SEXP a, SEXP i, SEXP j)
{
  int check;
  const int n = INTEGER(n_)[0];
  SEXP ret;
  PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  check = cosine_sparse_coo(INEDEX_FROM_1, n, LENGTH(a), REAL(a), INTEGER(i), INTEGER(j), REAL(ret));
  if (check)
    error("unable to allocate necessary memory");
  
  UNPROTECT(1);
  
  return ret;
}



SEXP R_sparsity_int(SEXP x)
{
  int i, j, count = 0;
  const int m = nrows(x), n = ncols(x);
  const double eps = 1.1 * DBL_EPSILON * DBL_MIN;
  int *x_pt = INTEGER(x);
  SEXP ret;
  
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      if (x_pt[i + m*j] == 0)
        count++;
    }
  }
  
  
  PROTECT(ret = allocVector(INTSXP, 1));
  INTEGER(ret)[0] = count;
  UNPROTECT(1);
  
  return ret;
}



SEXP R_sparsity_dbl(SEXP x, SEXP tol)
{
  int i, j, count = 0;
  const int m = nrows(x), n = ncols(x);
  const double eps = REAL(tol)[0];
  double *x_pt = REAL(x);
  SEXP ret;
  
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      if (fabs(x_pt[i + m*j]) < eps)
        count++;
    }
  }
  
  
  PROTECT(ret = allocVector(INTSXP, 1));
  INTEGER(ret)[0] = count;
  UNPROTECT(1);
  
  return ret;
}

