/*  Copyright (c) 2015-2016, Schmidt
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
#include "fastco.h"


#define CO_SIM 1
#define CO_ORR 2
#define CO_VAR 3

#define BADTYPE() error("Invalid 'type' argument; please report this to the package author")

// ---------------------------------------------
//  Dense
// ---------------------------------------------

SEXP R_co_mat(SEXP x, SEXP type_)
{
  SEXP ret;
  int check;
  const int type = INTEGER(type_)[0];
  const unsigned int m = nrows(x);
  const unsigned int n = ncols(x);
  PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  if (type == CO_SIM)
    check = cosine_mat(m, n, REAL(x), REAL(ret));
  else if (type == CO_ORR)
    check = pcor_mat(m, n, REAL(x), REAL(ret));
  else if (type == CO_VAR)
    check = covar_mat(m, n, REAL(x), REAL(ret));
  else
    BADTYPE();
  
  UNPROTECT(1);
  
  if (check)
    error("unable to allocate necessary memory");
  
  return ret;
}



SEXP R_co_vecvec(SEXP x, SEXP y, SEXP type_)
{
  SEXP ret;
  int check;
  const int type = INTEGER(type_)[0];
  const unsigned int n = LENGTH(x);
  double co;
  PROTECT(ret = allocVector(REALSXP, 1));
  
  if (type == CO_SIM)
    check = cosine_vecvec(n, REAL(x), REAL(y), &co);
  else if (type == CO_ORR)
    check = pcor_vecvec(n, REAL(x), REAL(y), &co);
  else if (type == CO_VAR)
    check = covar_vecvec(n, REAL(x), REAL(y), &co);
  else
    BADTYPE();
  
  REAL(ret)[0] = co;
  UNPROTECT(1);
  
  if (check)
    error("unable to allocate necessary memory");
  
  return ret;
}





// ---------------------------------------------
//  Sparse
// ---------------------------------------------

#define INEDEX_FROM_1 1

SEXP R_co_sparse(SEXP n_, SEXP a, SEXP i, SEXP j, SEXP type_)
{
  const int type = INTEGER(type_)[0];
  int check;
  const int n = INTEGER(n_)[0];
  SEXP ret;
  PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  if (type == CO_SIM)
    check = cosine_sparse_coo(INEDEX_FROM_1, n, LENGTH(a), REAL(a), INTEGER(i), INTEGER(j), REAL(ret));
  else
    BADTYPE();
  
  UNPROTECT(1);
  
  if (check)
    error("unable to allocate necessary memory");
  
  return ret;
}





// ---------------------------------------------
//  Sparse utils
// ---------------------------------------------

SEXP R_sparsity_int(SEXP x)
{
  int i, j, count = 0;
  const int m = nrows(x), n = ncols(x);
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
