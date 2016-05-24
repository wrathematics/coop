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

// Functions for computing covariance, (pearson) correlation, and cosine similarity

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "coop.h"
#include "omputils.h"
#include "lapack.h"


// ---------------------------------------------
//  Static utils
// ---------------------------------------------

// ddot replica using dgemm
static inline double ddot(const int n, const double * const restrict x, const double * const restrict y)
{
  const int one = 1;
  double dot;
  
  dgemm_(&(char){'t'}, &(char){'n'}, &one, &one, &n,
    &(double){1.0}, x, &n, y, &n, &(double){0.0}, &dot, &one);
  
  return dot;
}



// upper triangle of t(x) %*% x
static inline void crossprod(const int m, const int n, const double * const restrict x, double *restrict c)
{
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &(double){1.0}, x, &m, &(double){0.0}, c, &n);
}



// x[*, j] -= colmean(x[*, j])
static void remove_colmeans(const int m, const int n, double *restrict x)
{
  if (m == 0 || n == 0)
    return;
    
  const double div = 1. / ((double) m);
  
  #pragma omp parallel for default(none) shared(x) if(m*n > OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    double colmean = 0;
    
    // Get column mean
    SAFE_SIMD
    for (int i=0; i<m; i++)
      colmean += x[i   + m*j];
      
    colmean *= div;
    
    // Remove mean from column
    SAFE_SIMD
    for (int i=0; i<m; i++)
      x[i   + m*j] -= colmean;
  }
}



// compute the mean of a vector
static inline double mean(const int n, const double * const restrict x)
{
  const double divbyn = 1. / ((double) n);
  double mean = 0.;
  
  SAFE_FOR_SIMD
  for (int i=0; i<n; i++)
    mean += x[i];
  
  return mean*divbyn;
}





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
 * @param m,n
 * The number of rows/columns of the input matrix x.
 * @param x
 * The input mxn matrix.
 * @param cos
 * The output nxn matrix.
*/
int coop_cosine_mat(const int m, const int n, const double * const restrict x, double *restrict cos)
{
  crossprod(m, n, x, cos);
  coop_fill(n, cos);
  coop_symmetrize(n, cos);
  
  return 0;
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
int coop_cosine_vecvec(const int n, const double * const restrict x, const double * const restrict y, double *cos)
{
  double normx, normy;
  const double cp = ddot(n, x, y);
  
  crossprod(n, 1, x, &normx);
  crossprod(n, 1, y, &normy);
  
  *cos = cp / sqrt(normx * normy);
  return 0;
}





// ---------------------------------------------
//  Correlation
// ---------------------------------------------

/**
 * @brief
 * Compute the pearson correlation matrix.
 *
 * @details
 * The implementation is dominated by a symmetric rank-k update
 * via the BLAS function dsyrk().
 *
 * @param m,n
 * The number of rows/columns of the input matrix x.
 * @param x
 * The input mxn matrix.
 * @param cor
 * The output nxn matrix.
*/
int coop_pcor_mat(const int m, const int n, const double * const restrict x, double *restrict cor)
{
  double *x_cp = malloc(m*n*sizeof(*x));
  CHECKMALLOC(x_cp);
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  remove_colmeans(m, n, x_cp);
  crossprod(m, n, x_cp, cor);
  coop_fill(n, cor);
  coop_symmetrize(n, cor);
  
  free(x_cp);
  return 0;
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
int coop_pcor_vecvec(const int n, const double * const  const restrict x, const double * const  const restrict y, double *restrict cor)
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
  
  crossprod(n, 1, x_minusmean, &normx);
  crossprod(n, 1, y_minusmean, &normy);
  
  free(x_minusmean);
  free(y_minusmean);
  
  *cor = cp / sqrt(normx * normy);
  return 0;
}





// ---------------------------------------------
//  Covariance
// ---------------------------------------------


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
int coop_covar_mat(const int m, const int n, const double * const restrict x, double *restrict cov)
{
  double alpha = 1. / ((double) (m-1));
  double *x_cp = malloc(m*n*sizeof(*x));
  CHECKMALLOC(x_cp);
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  remove_colmeans(m, n, x_cp);
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &alpha, x_cp, &m, &(double){0.0}, cov, &n);
  coop_symmetrize(n, cov);
  
  free(x_cp);
  
  return 0;
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
int coop_covar_vecvec(const int n, const double * const restrict x, const double * const restrict y, double *restrict cov)
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
  return 0;
}
