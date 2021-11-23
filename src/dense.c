/*  Copyright (c) 2015-2016, 2021 Drew Schmidt
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


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils/safeomp.h"

#include "coop.h"
#include "utils/fill.h"
#include "utils/inverse.h"
#include "utils/mmult.h"
#include "utils/scale.h"
#include "utils/sumstats.h"
#include "utils/xpose.h"


// ---------------------------------------------
//  Cosine
// ---------------------------------------------

/**
 * @brief
 * Compute the cosine similarity matrix of a matrix.  This is
 * all pair-wise vector cosine similarities of the columns.
 *
 * @details
 * The implementation is dominated by a symmetric rank-k update
 * via the BLAS function dsyrk().
 *
 * @param trans
 * Perform cosine(x) or cosine(t(x))
 * @param m,n
 * The number of rows/columns of the input matrix x.
 * @param x
 * The input mxn matrix.
 * @param cos
 * The output nxn matrix.
*/
int coop_cosine_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cos)
{
  int ncols;
  
  ncols = n;
  crossprod(m, n, 1.0, x, m, cos);
  
  int ret = cosim_fill(ncols, cos);
  CHECKRET(ret);
  
  if (inv)
  {
    ret = inv_sym_chol(ncols, cos);
    CHECKRET(ret);
  }
  
  symmetrize(ncols, cos);
  
  return COOP_OK;
}

int coop_tcosine_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cos)
{
  int ncols;
  
  ncols = m;
  tcrossprod(m, n, 1.0, x, m, cos);
  
  int ret = cosim_fill(ncols, cos);
  CHECKRET(ret);
  
  if (inv)
  {
    ret = inv_sym_chol(ncols, cos);
    CHECKRET(ret);
  }
  
  symmetrize(ncols, cos);
  
  return COOP_OK;
}



int coop_cosine_matmat(const bool inv, const int m, const int nx,
  const double *const restrict x, const int ny, const double *const restrict y,
  double *restrict cos)
{
  matmult(true, false, 1.0, m, nx, x, m, ny, y, cos);
  
  double *diagx = malloc(nx * sizeof(*diagx));
  double *diagy = malloc(ny * sizeof(*diagy));
  CHECKMALLOC(diagx);
  CHECKMALLOC(diagy);
  
  SAFE_FOR_SIMD
  for (int i=0; i<nx; i++)
    crossprod(m, 1, 1.0, x+m*i, m, diagx+i);
  SAFE_FOR_SIMD
  for (int i=0; i<ny; i++)
    crossprod(m, 1, 1.0, y+m*i, m, diagy+i);
  
  for (int j=0; j<ny; j++)
  {
    for (int i=0; i<nx; i++)
      cos[i + nx*j] /= sqrt(diagx[i] * diagy[j]);
  }
  
  free(diagx);
  free(diagy);
  
  if (inv && nx == ny)
  {
    int ret = inv_gen_lu(nx, cos);
    CHECKRET(ret);
  }
  
  return COOP_OK;
}


int coop_tcosine_matmat(const bool inv, const int mx, const int n,
  const double *const restrict x, const int my, const double *const restrict y,
  double *restrict cos)
{
  matmult(false, true, 1.0, mx, n, x, my, n, y, cos);
  
  double *diagx = malloc(mx * sizeof(*diagx));
  double *diagy = malloc(my * sizeof(*diagy));
  CHECKMALLOC(diagx);
  CHECKMALLOC(diagy);
  
  SAFE_FOR_SIMD
  for (int i=0; i<mx; i++)
    tcrossprod(1, n, 1.0, x+i, mx, diagx+i);
  SAFE_FOR_SIMD
  for (int i=0; i<my; i++)
    tcrossprod(1, n, 1.0, y+i, my, diagy+i);
  
  for (int j=0; j<my; j++)
  {
    for (int i=0; i<mx; i++)
      cos[i + mx*j] /= sqrt(diagx[i] * diagy[j]);
  }
  
  free(diagx);
  free(diagy);
  
  if (inv && mx == my)
  {
    int ret = inv_gen_lu(mx, cos);
    CHECKRET(ret);
  }
  
  return COOP_OK;
}



/**
 * @brief
 * Compute the cosine similarity between two vectors.
 *
 * @details
 * The implementation uses a dgemm() to compute the dot product
 * of x and y, and then two dsyrk() calls to compute the (square of)
 * the norms of x and y.
 *
 * @param n
 * The length of the x and y vectors.
 * @param x,y
 * The input vectors.
 *
 * @return
 * The cosine similarity between the two vectors.
*/
int coop_cosine_vecvec(const int n, const double *const restrict x,
  const double *const restrict y, double *cos)
{
  double normx, normy;
  const double cp = ddot(n, x, y);
  
  crossprod(n, 1, 1.0, x, n, &normx);
  crossprod(n, 1, 1.0, y, n, &normy);
  
  *cos = cp / sqrt(normx * normy);
  return COOP_OK;
}



// ---------------------------------------------
// Pearson Correlation
// ---------------------------------------------

static inline int coop_pcor_mat_work(const bool inv, const int m, const int n,
  double *const restrict x, double *restrict cor)
{
  remove_colmeans(m, n, x);
  crossprod(m, n, 1.0, x, m, cor);
  free(x);
  
  int ret = cosim_fill(n, cor);
  CHECKRET(ret);
  
  if (inv)
  {
    ret = inv_sym_chol(n, cor);
    CHECKRET(ret);
  }
  
  symmetrize(n, cor);
  
  return COOP_OK;
}

/**
 * @brief
 * Compute the pearson correlation matrix.
 *
 * @details
 * The implementation is dominated by a symmetric rank-k update
 * via the BLAS function dsyrk().
 *
 * @param inv
 * Invert after computing?
 * @param m,n
 * The number of rows/columns of the input matrix x.
 * @param x
 * The input mxn matrix.
 * @param cor
 * The output correlation matrix.
*/
int coop_pcor_mat(const bool inv, const int m, const int n,
  const double * const restrict x, double *restrict cor)
{
  double *x_cp = malloc(m*n*sizeof(*x));
  CHECKMALLOC(x_cp);
  
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  coop_pcor_mat_work(inv, m, n, x_cp, cor);
  
  return COOP_OK;
}

int coop_tpcor_mat(const bool inv, const int m, const int n,
  const double * const restrict x, double *restrict cor)
{
  double *x_cp = malloc(m*n*sizeof(*x));
  CHECKMALLOC(x_cp);
  
  xpose(m, n, x, x_cp);
  
  coop_pcor_mat_work(inv, n, m, x_cp, cor);
  
  return COOP_OK;
}



/**
 * pcor(x, y)
 * 
 * @brief
 * Compute the pearson correlation matrix.
 *
 * @details
 * The implementation is dominated by a symmetric rank-k update
 * via the BLAS function dsyrk().
 *
 * @param trans
 * Transpose before computing?
 * @param inv
 * Invert after computing?
 * @param m,n
 * The number of rows/columns of the input matrix x.
 * @param x,y
 * The input matrices.
 * @param cor
 * The output correlation matrix.
*/
int coop_pcor_matmat(const bool inv, const int m, const int nx,
  const double *const restrict x, const int ny, const double *const restrict y,
  double *restrict cor)
{
  int ret = 0;
  double *x_cp = malloc(m*nx * sizeof(*x));
  double *y_cp = malloc(m*ny * sizeof(*y));
  CHECKMALLOC(x_cp);
  CHECKMALLOC(y_cp);
  
  memcpy(x_cp, x, m*nx*sizeof(*x));
  memcpy(y_cp, y, m*ny*sizeof(*y));
  
  scale_nostore(true, true, m, nx, x_cp);
  scale_nostore(true, true, m, ny, y_cp);
  
  const double alpha = 1. / ((double) (m-1));
  
  matmult(true, false, alpha, m, nx, x_cp, m, ny, y_cp, cor);
  free(x_cp);
  free(y_cp);
  
  if (inv)
    ret = inv_sym_chol(m, cor);
  
  return ret;
}



int coop_tpcor_matmat(const bool inv, const int mx, const int n,
  const double *const restrict x, const int my, const double *const restrict y,
  double *restrict cor)
{
  int ret = 0;
  double *x_cp = malloc(mx*n * sizeof(*x));
  double *y_cp = malloc(my*n * sizeof(*y));
  CHECKMALLOC(x_cp);
  CHECKMALLOC(y_cp);
  
  xpose(mx, n, x, x_cp);
  xpose(my, n, y, y_cp);
  
  scale_nostore(true, true, n, mx, x_cp);
  scale_nostore(true, true, n, my, y_cp);
  
  const double alpha = 1. / ((double) (n-1));
  
  matmult(true, false, alpha, n, mx, x_cp, n, my, y_cp, cor);
  free(x_cp);
  free(y_cp);
  
  if (inv)
    ret = inv_sym_chol(n, cor);
  
  return ret;
}



/**
 * @brief
 * Compute the pearson correlation between two vectors.
 *
 * @details
 * The implementation uses a dgemm() to compute the dot product
 * of x and y, and then two dsyrk() calls to compute the (square of)
 * the norms of x and y.
 *
 * @param n
 * The length of the x and y vectors.
 * @param x,y
 * The input vectors.
 *
 * @return
 * The correlation between the two vectors.
*/
int coop_pcor_vecvec(const int n, const double *const restrict x,
  const double *const restrict y, double *restrict cor)
{
  double normx, normy;
  
  double *x_minusmean = malloc(n*sizeof(*x));
  CHECKMALLOC(x_minusmean);
  double *y_minusmean = malloc(n*sizeof(*y));
  CHECKMALLOC(y_minusmean);
  
  const double meanx = mean(n, x);
  const double meany = mean(n, y);
  
  SAFE_PARALLEL_FOR_SIMD
  for (int i=0; i<n; i++)
  {
    x_minusmean[i] = x[i] - meanx;
    y_minusmean[i] = y[i] - meany;
  }
  
  const double cp = ddot(n, x_minusmean, y_minusmean);
  
  crossprod(n, 1, 1.0, x_minusmean, n, &normx);
  crossprod(n, 1, 1.0, y_minusmean, n, &normy);
  
  free(x_minusmean);
  free(y_minusmean);
  
  *cor = cp / sqrt(normx * normy);
  return COOP_OK;
}



// ---------------------------------------------
//  Covariance
// ---------------------------------------------

static inline int coop_covar_mat_work(const bool inv, const int m, const int n,
  double *const restrict x, double *restrict cov)
{
  const double alpha = 1. / ((double) (m-1));
  
  remove_colmeans(m, n, x);
  crossprod(m, n, alpha, x, m, cov);
  free(x);
  
  if (inv)
  {
    int ret = inv_sym_chol(n, cov);
    CHECKRET(ret);
  }
  
  symmetrize(n, cov);
  
  return COOP_OK;
}

/**
 * @file
 * @brief Covariance.
 *
 * @details
 * Computes the variance-covariance matrix.  Centering is done in-place.
 *
 * @param method
 * Input.  The form the covariance matrix takes (pearson, kendall,
 * spearman).  Currently only pearson works.
 * @param m,n
 * Inputs.  Problem size (dims of x)
 * @param x
 * Input.  The data matrix.
 * @param coc
 * Output.  The covariance matrix.
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int coop_covar_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cov)
{
  double *x_cp = malloc(m*n*sizeof(*x));
  CHECKMALLOC(x_cp);
  
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  return coop_covar_mat_work(inv, m, n, x_cp, cov);
}

int coop_tcovar_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cov)
{
  double *x_cp = malloc(m*n*sizeof(*x));
  CHECKMALLOC(x_cp);
  
  xpose(m, n, x, x_cp);
  
  return coop_covar_mat_work(inv, n, m, x_cp, cov);
}



// covar(x,y)
int coop_covar_matmat(const bool trans, const bool inv, const int m,
  const int n, const double *const restrict x, 
  const double *const restrict y, double *restrict cov)
{
  int ret = 0;
  int nrows, ncols;
  double *x_cp = malloc(m*n * sizeof(*x));
  CHECKMALLOC(x_cp);
  double *y_cp = malloc(m*n * sizeof(*y));
  if (y_cp == NULL)
  {
    free(x_cp);
    return COOP_BADMALLOC;
  }
  
  
  if (trans)
  {
    xpose(m, n, x, x_cp);
    xpose(m, n, y, y_cp);
    nrows = n;
    ncols = m;
  }
  else
  {
    memcpy(x_cp, x, m*n*sizeof(*x));
    memcpy(y_cp, y, m*n*sizeof(*y));
    nrows = m;
    ncols = n;
  }
  
  
  const double alpha = 1. / ((double) (nrows-1));
  
  //TODO FIXME make tremove_colmeans and use the BLAS more efficiently...
  remove_colmeans(nrows, ncols, x_cp);
  remove_colmeans(nrows, ncols, y_cp);
  
  // matmult(!trans, trans, alpha, nrows, ncols, x_cp, nrows, ncols, y_cp, cov);
  matmult(true, false, alpha, nrows, ncols, x_cp, nrows, ncols, y_cp, cov);
  free(x_cp);
  free(y_cp);
  
  if (inv)
    ret = inv_sym_chol(ncols, cov);
  
  return ret;
}



/**
 * @brief
 * Compute the covariance between two vectors.
 *
 * @details
 * The implementation uses a dgemm() to compute the dot product
 * of x and y, and then two dsyrk() calls to compute the (square of)
 * the norms of x and y.
 *
 * @param n
 * The length of the x and y vectors.
 * @param x,y
 * The input vectors.
 *
 * @return
 * The variance of the vectors.
*/
int coop_covar_vecvec(const int n, const double *const restrict x,
  const double *const restrict y, double *restrict cov)
{
  const double recip_n = (double) 1. / (n-1);
  double sum_xy = 0., sum_x = 0., sum_y = 0.;
  
  #ifdef OMP_VER_4
  #pragma omp simd reduction(+: sum_xy, sum_x, sum_y)
  #endif
  for (int i=0; i<n; i++)
  {
    const double tx = x[i];
    const double ty = y[i];
    
    sum_xy += tx*ty;
    sum_x += tx;
    sum_y += ty;
  }
  
  *cov = (sum_xy - (sum_x*sum_y*((double) 1./n))) * recip_n;
  return COOP_OK;
}
