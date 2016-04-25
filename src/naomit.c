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

#include "omputils.h"

#define INT(x) INTEGER(x)

#define COPYVEC(in,out,len,TYPE) \
  PROTECT(out = allocVector(TYPE##SXP, len)); \
  memcpy(TYPE(out), TYPE(in), len*sizeof(TYPE(in)));

#define COPYMAT(in,out,m,n,TYPE) \
  PROTECT(out = allocMatrix(TYPE##SXP, m, n)); \
  memcpy(TYPE(out), in, m*n*sizeof(in)); \
  UNPROTECT(1);


// --------------------------------------------------------------
// dense
// --------------------------------------------------------------

// faster to index each element and operate accordingly, but
// this is too memory expensive for most applications
// note: R does this anyway because, well, R...
static SEXP R_fast_naomit_dbl_small(const int m, const int n, const double *x)
{
  SEXP ret;
  const int len = m*n;
  int i, j, mj;
  int *na_vec_ind = (int*) calloc(len, sizeof(*na_vec_ind));
  int m_fin = m;
  int row;

  // get indices of NA's
  SAFE_FOR_SIMD
  for (i=0; i<len; i++)
  {
    if (ISNA(x[i]) || ISNAN(x[i]))
      na_vec_ind[i] = 1;
  }

  // adjust col index; turn first column of the NA indices
  // to track which rows should go
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
    COPYMAT(x, ret, m, n, REAL);
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



static SEXP R_fast_naomit_dbl_big(const int m, const int n, const double *x)
{
  SEXP ret;
  int i, j, mj;
  int *rows = (int*) calloc(m, sizeof(*rows));
  int m_fin = m;
  int row;

  // get indices of NA's
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

  // get number of rows of output
  SAFE_FOR_SIMD
  for (i=0; i<m; i++)
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



static SEXP R_fast_naomit_int_small(const int m, const int n, const int *x)
{
  SEXP ret;
  const int len = m*n;
  int i, j, mj;
  int *na_vec_ind = (int*) calloc(len, sizeof(*na_vec_ind));
  int m_fin = m;
  int row;

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
    COPYMAT(x, ret, m, n, INT);
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



static SEXP R_fast_naomit_int_big(const int m, const int n, const int *x)
{
  SEXP ret;
  int i, j, mj;
  int *rows = (int*) calloc(m, sizeof(*rows));
  int m_fin = m;
  int row;

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
    COPYMAT(x, ret, m, n, INT);
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
  if (isReal(x))
    return R_fast_naomit_dbl(x);
  else if (isInteger(x))
    return R_fast_naomit_int(x);
  else
    error("'x' must be numeric");
}



SEXP R_naomit_vecvec(SEXP x_, SEXP y_)
{
  int i;
  const int n = LENGTH(x_);
  double *x, *y;
  SEXP x_ret, y_ret;
  SEXP ret;

  // COPYVEC(x_, x_ret, n, REAL);
  // COPYVEC(y_, y_ret, n, REAL);
  // double *x = REAL(x_);
  // double *y = REAL(y_);

  x = malloc(n * sizeof(*x));
  memcpy(x, REAL(x_), n*sizeof(*x));
  y = malloc(n * sizeof(*y));
  memcpy(y, REAL(y_), n*sizeof(*y));

  for (i=0; i<n; i++)
  {
    if (ISNA(x[i]) || ISNAN(x[i]))
      y[i] = x[i];
    else if (ISNA(y[i]) || ISNAN(y[i]))
      x[i] = y[i];
  }

  x_ret = R_fast_naomit_dbl_small(n, 1, x);
  y_ret = R_fast_naomit_dbl_small(n, 1, y);

  free(x);
  free(y);

  PROTECT(ret = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ret, 0, x_ret);
  SET_VECTOR_ELT(ret, 1, y_ret);

  UNPROTECT(1);
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
