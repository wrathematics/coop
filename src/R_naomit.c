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

// NOTE valgrind reports "possibly lost" memory errors; this is GNU OMP's fault:
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=36298


#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

#include "utils/safeomp.h"

#define INT(x) INTEGER(x)

#define COPYVEC(in,out,len,TYPE) \
  PROTECT(out = allocVector(TYPE##SXP, len)); \
  memcpy(TYPE(out), TYPE(in), len*sizeof(*TYPE(in))); 

#define COPYMAT(in,out,m,n,TYPE) \
  PROTECT(out = allocMatrix(TYPE##SXP, m, n)); \
  memcpy(TYPE(out), in, m*n*sizeof(*in)); \
  UNPROTECT(1);

#define THROW_MEMERR error("unable to allocate necessary memory")
#define R_CHECKMALLOC(ptr) if (ptr == NULL) THROW_MEMERR


// --------------------------------------------------------------
// dense
// --------------------------------------------------------------

// faster to index each element and operate accordingly, but
// this is too memory expensive for most applications
// note: R does this anyway because, well, R...
static SEXP R_fast_naomit_dbl_small(const int m, const int n, const double *const x)
{
  SEXP ret;
  const int len = m*n;
  int m_fin = m;
  int *na_vec_ind = calloc(len, sizeof(*na_vec_ind));
  R_CHECKMALLOC(na_vec_ind);
  
  
  // get indices of NA's
  PLEASE_VECTORIZE
  for (int i=0; i<len; i++)
  {
    if (ISNA(x[i]) || ISNAN(x[i]))
      na_vec_ind[i] = 1;
  }
  
  // adjust col index; turn first column of the NA indices
  // to track which rows should go
  for (int j=1; j<n; j++)
  {
    const int mj = m*j;
    
    for (int i=0; i<m; i++)
    {
      if (na_vec_ind[i + mj])
        na_vec_ind[i] = 1;
    }
  }
  
  // get number of rows of output
  for (int i=0; i<m; i++)
    m_fin -= na_vec_ind[i];
  
  // do a cheap copy if the matrix is identical
  if (m_fin == m)
  {
    COPYMAT(x, ret, m, n, REAL);
    free(na_vec_ind);
    return ret;
  }
  
  // build reduced matrix
  PROTECT(ret = allocMatrix(REALSXP, m_fin, n));
  double *retptr = REAL(ret);
  
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    int row = 0;
    
    for (int i=0; i<m; i++)
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



static SEXP R_fast_naomit_dbl_big(const int m, const int n, const double *const x)
{
  SEXP ret;
  int m_fin = m;
  int *rows = calloc(m, sizeof(*rows));
  R_CHECKMALLOC(rows);
  
  // get indices of NA's
  #pragma omp parallel for shared(rows)
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    
    for (int i=0; i<m; i++)
    {
      if (ISNAN(x[i + mj])) // ISNAN for doubles will return true for NA and NaN
        rows[i] = 1;
    }
  }
  
  // get number of rows of output
  for (int i=0; i<m; i++)
    m_fin -= rows[i];
  
  // do a cheap copy if the matrix is identical
  if (m_fin == m)
  {
    COPYMAT(x, ret, m, n, REAL);
    free(rows);
    return ret;
  }
  
  PROTECT(ret = allocMatrix(REALSXP, m_fin, n));
  double *retptr = REAL(ret);
  
  // build reduced matrix
  #pragma omp parallel for shared(rows, retptr, m_fin)
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    int row = 0;
    
    for (int i=0; i<m; i++)
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



SEXP R_fast_naomit_dbl(SEXP x_)
{
  const int m = nrows(x_);
  const int n = ncols(x_);
  const double *x = REAL(x_);
  
  if (m*n < OMP_MIN_SIZE)
    return R_fast_naomit_dbl_small(m, n, x);
  else
    return R_fast_naomit_dbl_big(m, n, x);
}



static SEXP R_fast_naomit_int_small(const int m, const int n, const int *const x)
{
  SEXP ret;
  const int len = m*n;
  int m_fin = m;
  int *na_vec_ind = calloc(len, sizeof(*na_vec_ind));
  R_CHECKMALLOC(na_vec_ind);
  
  // get indices of NA's
  PLEASE_VECTORIZE
  for (int i=0; i<len; i++)
  {
    if (x[i] == NA_INTEGER)
      na_vec_ind[i] = 1;
  }
  
  // adjust col index
  for (int j=1; j<n; j++)
  {
    const int mj = m*j;
    
    for (int i=0; i<m; i++)
    {
      if (na_vec_ind[i + mj])
        na_vec_ind[i] = 1;
    }
  }
  
  // get number of rows of output
  for (int i=0; i<m; i++)
    m_fin -= na_vec_ind[i];
  
  // do a cheap copy if the matrix is identical
  if (m_fin == m)
  {
    COPYMAT(x, ret, m, n, INT);
    free(na_vec_ind);
    return ret;
  }
  
  // build reduced matrix
  PROTECT(ret = allocMatrix(INTSXP, m_fin, n));
  int *retptr = INTEGER(ret);
  
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    int row = 0;
    
    for (int i=0; i<m; i++)
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



static SEXP R_fast_naomit_int_big(const int m, const int n, const int *const x)
{
  SEXP ret;
  int m_fin = m;
  int *rows = calloc(m, sizeof(*rows));
  R_CHECKMALLOC(rows);
  
  #pragma omp parallel for shared(rows, NA_INTEGER)
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    
    for (int i=0; i<m; i++)
    {
      if (x[i + mj] == NA_INTEGER)
        rows[i] = 1;
    }
  }
  
  for (int i=0; i<m; i++)
    m_fin -= rows[i];
  
  if (m_fin == m)
  {
    COPYMAT(x, ret, m, n, INT);
    free(rows);
    return ret;
  }
  
  PROTECT(ret = allocMatrix(INTSXP, m_fin, n));
  int *retptr = INTEGER(ret);
  
  #pragma omp parallel for shared(rows, retptr, m_fin)
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    int row = 0;
    
    for (int i=0; i<m; i++)
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



SEXP R_fast_naomit_int(SEXP x_)
{
  const int m = nrows(x_);
  const int n = ncols(x_);
  
  const int *x = INTEGER(x_);
  
  if (m*n < OMP_MIN_SIZE)
    return R_fast_naomit_int_small(m, n, x);
  else
    return R_fast_naomit_int_big(m, n, x);
}



SEXP R_fast_naomit(SEXP x)
{
  switch (TYPEOF(x))
  {
    case REALSXP:
      return R_fast_naomit_dbl(x);
    
    case INTSXP:
      return R_fast_naomit_int(x);
    
    default:
      error("'x' must be numeric");
  }
}



SEXP R_naomit_vecvec(SEXP x_, SEXP y_)
{
  const int n = LENGTH(x_);
  double *x, *y;
  SEXP x_ret, y_ret;
  SEXP ret;
  
  // COPYVEC(x_, x_ret, n, REAL);
  // COPYVEC(y_, y_ret, n, REAL);
  // double *x = REAL(x_);
  // double *y = REAL(y_);
  
  x = malloc(n * sizeof(*x));
  R_CHECKMALLOC(x);
  
  y = malloc(n * sizeof(*y));
  if (y == NULL)
  {
    free(x);
    THROW_MEMERR;
  }
  
  memcpy(x, REAL(x_), n*sizeof(*x));
  memcpy(y, REAL(y_), n*sizeof(*y));
  
  for (int i=0; i<n; i++)
  {
    if (ISNA(x[i]) || ISNAN(x[i]))
      y[i] = x[i];
    else if (ISNA(y[i]) || ISNAN(y[i]))
      x[i] = y[i];
  }
  
  PROTECT(x_ret = R_fast_naomit_dbl_small(n, 1, x));
  PROTECT(y_ret = R_fast_naomit_dbl_small(n, 1, y));
  
  free(x);
  free(y);
  
  PROTECT(ret = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ret, 0, x_ret);
  SET_VECTOR_ELT(ret, 1, y_ret);
  
  UNPROTECT(3);
  return ret;
}



// --------------------------------------------------------------
// sparse
// --------------------------------------------------------------

/*
SEXP R_naomit_coo(SEXP a_in_, SEXP i_in_, SEXP j_in_)
{
  SEXP ret;
  SEXP a_out_, i_out_, j_out_;
  const double *a_in = REAL(a_in_);
  const int *i_in = INTEGER(i_in_);
  const int *j_in = INTEGER(j_in_);
  const int len_in = LENGTH(a_in_);
  int k;
  int len_out = 0;
  
  
  // Find all NA rows
  for (k=0; k<len_in; k++)
  {
  
  }
  
  
  
  
  // build reduced matrix
  PROTECT(a_out_ = allocVector(REALSXP, len_out));
  PROTECT(i_out_ = allocVector(INTSXP, len_out));
  PROTECT(j_out_ = allocVector(INTSXP, len_out));
  int *i_out = INTEGRE(i_out_);
  int *j_out = INTEGRE(j_out_);
  
  
  
  
  
  // Set return
  PROTECT(ret = allocVector(VECSXP, 3));
  
  SET_VECTOR_ELT(ret, 0, a_out_);
  SET_VECTOR_ELT(ret, 1, i_out_);
  SET_VECTOR_ELT(ret, 2, j_out_);
  
  UNPROTECT(4);
  return ret;
}
*/
