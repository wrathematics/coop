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

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

#include "omp.h"


static SEXP copymat_dbl(const int m, const int n, SEXP x_)
{
  SEXP ret;
  const double *x = REAL(x_);
  PROTECT(ret = allocMatrix(REALSXP, m, n));
  double *retptr = REAL(ret);
  memcpy(retptr, x, m*n*sizeof(*x));
  
  UNPROTECT(1);
  return ret;
}



static SEXP copymat_int(const int m, const int n, SEXP x_)
{
  SEXP ret;
  const int *x = INTEGER(x_);
  PROTECT(ret = allocMatrix(INTSXP, m, n));
  int *retptr = INTEGER(ret);
  memcpy(retptr, x, m*n*sizeof(*x));
  
  UNPROTECT(1);
  return ret;
}



// faster to index each element and operate accordingly, but
// this is too memory expensive for most applications
// note: R does this anyway because, well, R...
static SEXP R_fast_naomit_dbl_small(const int m, const int n, SEXP x_)
{
  SEXP ret;
  const int len = m*n;
  int i, j, mj;
  int *na_vec_ind = (int*) calloc(len, sizeof(*na_vec_ind));
  int m_fin = m;
  int row;
  
  const double *x = REAL(x_);
  
  // get indices of NA's
  SAFE_FOR_SIMD
  for (i=0; i<len; i++)
  {
    if (ISNA(x[i]) || ISNAN(x[i]))
      na_vec_ind[i] = 1;
  }
  
  // adjust col index
  for (j=1; j<n; j++)
  {
    mj = m*j;
    for (i=0; i<m; i++)
    {
      if (na_vec_ind[i + mj])
      {
        na_vec_ind[i] = 1;
      }
    }
  }
  
  // get number of rows of output
  SAFE_FOR_SIMD
  for (i=0; i<m; i++)
    m_fin -= na_vec_ind[i];
  
  // do a cheap copy if the matrix is identical
  if (m_fin == m)
  {
    ret = copymat_dbl(m, n, x_);
    free(na_vec_ind);
    return ret;
  }
  
  // build reduced matrix
  PROTECT(ret = allocMatrix(REALSXP, m_fin, n));
  double *retptr = REAL(ret);
  
  SAFE_FOR_SIMD
  for (j=0; j<n; j++)
  {
    mj = m*j;
    row = 0;
    
    for (i=0; i<m; i++)
    {
      if (!na_vec_ind[i%m])
      {
        retptr[row + m_fin*j] = x[i + mj];
        row++;
      }
    }
  }
  
  free(na_vec_ind);
  UNPROTECT(1);
  return ret;
}



static SEXP R_fast_naomit_dbl_big(const int m, const int n, SEXP x_)
{
  SEXP ret;
  int i, j, mj;
  int *rows = (int*) calloc(m, sizeof(*rows));
  int m_fin = m;
  int row;
  
  const double *x = REAL(x_);
  
  #pragma omp parallel for default(shared) private(i, j, mj)
  for (j=0; j<n; j++)
  {
    mj = m*j;
    
    SAFE_SIMD
    for (i=0; i<m; i++)
    {
      if (ISNA(x[i + mj]) || ISNAN(x[i + mj]))
        rows[i] = 1;
    }
  }
  
  SAFE_FOR_SIMD
  for (i=0; i<m; i++)
    m_fin -= rows[i];
  
  if (m_fin == m)
  {
    ret = copymat_dbl(m, n, x_);
    free(rows);
    return ret;
  }
  
  PROTECT(ret = allocMatrix(REALSXP, m_fin, n));
  double *retptr = REAL(ret);
  
  #pragma omp parallel for default(shared) private(i, j, row, mj)
  for (j=0; j<n; j++)
  {
    mj = m*j;
    row = 0;
    
    SAFE_SIMD
    for (i=0; i<m; i++)
    {
      if (!rows[i])
      {
        retptr[row + m_fin*j] = x[i + mj];
        row++;
      }
    }
  }
  
  free(rows);
  UNPROTECT(1);
  return ret;
}



SEXP R_fast_naomit_dbl(SEXP x)
{
  const int m = nrows(x);
  const int n = ncols(x);
  
  if (m*n < OMP_MIN_SIZE)
    return R_fast_naomit_dbl_small(m, n, x);
  else
    return R_fast_naomit_dbl_big(m, n, x);
}



static SEXP R_fast_naomit_int_small(const int m, const int n, SEXP x_)
{
  SEXP ret;
  const int len = m*n;
  int i, j, mj;
  int *na_vec_ind = (int*) calloc(len, sizeof(*na_vec_ind));
  int m_fin = m;
  int row;
  
  const int *x = INTEGER(x_);
  
  // get indices of NA's
  SAFE_FOR_SIMD
  for (i=0; i<len; i++)
  {
    if (x[i] == NA_INTEGER)
      na_vec_ind[i] = 1;
  }
  
  // adjust col index
  for (j=1; j<n; j++)
  {
    mj = m*j;
    for (i=0; i<m; i++)
    {
      if (na_vec_ind[i + mj])
      {
        na_vec_ind[i] = 1;
      }
    }
  }
  
  // get number of rows of output
  SAFE_FOR_SIMD
  for (i=0; i<m; i++)
    m_fin -= na_vec_ind[i];
  
  // do a cheap copy if the matrix is identical
  if (m_fin == m)
  {
    ret = copymat_int(m, n, x_);
    free(na_vec_ind);
    return ret;
  }
  
  // build reduced matrix
  PROTECT(ret = allocMatrix(INTSXP, m_fin, n));
  int *retptr = INTEGER(ret);
  
  SAFE_FOR_SIMD
  for (j=0; j<n; j++)
  {
    mj = m*j;
    row = 0;
    
    for (i=0; i<m; i++)
    {
      if (!na_vec_ind[i%m])
      {
        retptr[row + m_fin*j] = x[i + mj];
        row++;
      }
    }
  }
  
  free(na_vec_ind);
  UNPROTECT(1);
  return ret;
}



static SEXP R_fast_naomit_int_big(const int m, const int n, SEXP x_)
{
  SEXP ret;
  int i, j, mj;
  int *rows = (int*) calloc(m, sizeof(*rows));
  int m_fin = m;
  int row;
  
  const int *x = INTEGER(x_);
  
  #pragma omp parallel for default(shared) private(i, j, mj)
  for (j=0; j<n; j++)
  {
    mj = m*j;
    
    SAFE_SIMD
    for (i=0; i<m; i++)
    {
      if (x[i + mj] == NA_INTEGER)
        rows[i] = 1;
    }
  }
  
  SAFE_FOR_SIMD
  for (i=0; i<m; i++)
    m_fin -= rows[i];
  
  if (m_fin == m)
  {
    ret = copymat_int(m, n, x_);
    free(rows);
    return ret;
  }
  
  PROTECT(ret = allocMatrix(INTSXP, m_fin, n));
  int *retptr = INTEGER(ret);
  
  #pragma omp parallel for default(shared) private(i, j, row, mj)
  for (j=0; j<n; j++)
  {
    mj = m*j;
    row = 0;
    
    SAFE_SIMD
    for (i=0; i<m; i++)
    {
      if (!rows[i])
      {
        retptr[row + m_fin*j] = x[i + mj];
        row++;
      }
    }
  }
  
  free(rows);
  UNPROTECT(1);
  return ret;
}



SEXP R_fast_naomit_int(SEXP x)
{
  const int m = nrows(x);
  const int n = ncols(x);
  
  if (m*n < OMP_MIN_SIZE)
    return R_fast_naomit_int_small(m, n, x);
  else
    return R_fast_naomit_int_big(m, n, x);
}



SEXP R_fast_naomit(SEXP x)
{
  if (isReal(x))
    return R_fast_naomit_dbl(x);
  else if (isInteger(x))
    return R_fast_naomit_int(x);
  else
    error("'x' must be numeric");
}
