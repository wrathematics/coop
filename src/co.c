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


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cosine.h"
#include "omp.h"


// ---------------------------------------------
//  Dense
// ---------------------------------------------

// ddot replica using dgemm
static inline double ddot(const int n, const double *restrict x, const double *restrict y)
{
  const int one = 1;
  double dot;
  
  dgemm_(&(char){'t'}, &(char){'n'}, &one, &one, &n, 
    &(double){1.0}, x, &n, y, &n, &(double){0.0}, &dot, &one);
  
  return dot;
}



// upper triangle of t(x) %*% x
static inline void crossprod(const int m, const int n, const double *restrict x, double *restrict c)
{
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &(double){1.0}, x, &m, &(double){0.0}, c, &n);
}



static inline void diag2one(const unsigned int n, double *restrict x)
{
  int i;
  
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
    x[i + n*i] = 1.0;
}



// replaces upper triangle of the crossproduct of a matrix with its cosine similarity
static inline void cosim_fill(const unsigned int n, double *restrict cp)
{
  int i, j;
  double diagj;
  
  #pragma omp parallel for private(i,j,diagj) default(shared) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    diagj = cp[j + n*j];
    
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      cp[i + n*j] /= sqrt(cp[i + n*i] * diagj);
  }
  
  diag2one(n, cp);
}



// Copy lower triangle to upper
static inline void symmetrize(const int n, double *restrict x)
{
  int i, j;
  
  #pragma omp parallel for private(i) default(shared) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      x[j + n*i] = x[i + n*j];
  }
}



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
void cosine_mat(const int m, const int n, const double *restrict x, double *restrict cos)
{
  crossprod(m, n, x, cos);
  cosim_fill(n, cos);
  symmetrize(n, cos);
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
double cosine_vecvec(const int n, const double *restrict x, const double *restrict y)
{
  double normx, normy;
  const double cp = ddot(n, x, y);
  
  crossprod(n, 1, x, &normx);
  crossprod(n, 1, y, &normy);
  
  return cp / sqrt(normx * normy);
}



static void remove_colmeans(const int m, const int n, double *restrict x)
{
  int i, j;
  double colmean;
  
  if (m == 0 || n == 0) 
    return;
  
  const double div = 1. / ((double) m);
  
  #pragma omp parallel for private(i, j, colmean) shared(x) if(m*n > OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    colmean = 0;
    
    // Get column mean
    SAFE_SIMD
    for (i=0; i<m; i++)
      colmean += x[i   + m*j];
    
    colmean *= div;
    
    // Remove mean from column
    SAFE_SIMD
    for (i=0; i<m; i++)
      x[i   + m*j] -= colmean;
    
  }
}



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
void pcor_mat(const int m, const int n, const double *restrict x, double *restrict cor)
{
  double *x_cp = malloc(m*n*sizeof(*x));
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  remove_colmeans(m, n, x_cp);
  crossprod(m, n, x_cp, cor);
  cosim_fill(n, cor);
  symmetrize(n, cor);
  
  free(x_cp);
}



static inline double mean(const int n, const double *restrict x)
{
  int i;
  const double divbyn = 1. / ((double) n);
  double mean = 0.;
  
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
    mean += x[i];
  
  return mean*divbyn;
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
double pcor_vecvec(const int n, const double *restrict x, const double *restrict y)
{
  int i;
  double normx, normy;
  double *x_minusmean = malloc(n*sizeof(*x));
  double *y_minusmean = malloc(n*sizeof(*y));
  
  const double meanx = mean(n, x);
  const double meany = mean(n, y);
  
  SAFE_PARALLEL_FOR_SIMD
  for (i=0; i<n; i++)
  {
    x_minusmean[i] = x[i] - meanx;
    y_minusmean[i] = y[i] - meany;
  }
  
  const double cp = ddot(n, x_minusmean, y_minusmean);
  
  crossprod(n, 1, x_minusmean, &normx);
  crossprod(n, 1, y_minusmean, &normy);
  
  free(x_minusmean);
  free(y_minusmean);
  
  return cp / sqrt(normx * normy);
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
void covar_mat(const int m, const int n, const double *restrict x, double *restrict cov)
{
  int info = 0;
  double alpha = 1. / ((double) (m-1));
  double *x_cp = malloc(m*n*sizeof(*x));
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  remove_colmeans(m, n, x_cp);
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &alpha, x_cp, &m, &(double){0.0}, cov, &n);
  symmetrize(n, cov);
  
  free(x_cp);
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
double covar_vecvec(const int n, const double *restrict x, const double *restrict y)
{
  const double recip_n = (double) 1. / (n-1);
  double sum_xy = 0., sum_x = 0., sum_y = 0.;
  double tx, ty;
  
  #pragma omp simd reduction(+: sum_xy, sum_x, sum_y)
  for (int i=0; i<n; i++)
  {
    tx = x[i];
    ty = y[i];
    
    sum_xy += tx*ty;
    sum_x += tx;
    sum_y += ty;
  }
  
  return (sum_xy - (sum_x*sum_y*((double) 1./n))) * recip_n;
}





// ---------------------------------------------
//  Sparse
// ---------------------------------------------

static inline void set2zero(const unsigned int n, double *restrict x)
{
  memset(x, 0.0, n*sizeof(*x));
}



// NaN-out a row/column of cos matrix for numerical compatibility with dense methods
static inline void set2nan(const int j, const int n, double *restrict cos)
{
  int i;
  
  for (i=j; i<n; i++)
    cos[i + n*j] = NAN;
  
  for (i=0; i<j; i++)
    cos[j + n*i] = NAN;
}



static inline double sparsedot_self(const int vecstart, const int vecend, const int *rows, const double *a)
{
  int i;
  double dot = 0.0;
  
  for (i=vecstart; i<=vecend; i++)
    dot += a[i]*a[i];
  
  return dot;
}



// get the first and last indices in the COO for column i
static inline void get_startend(const int len, const int ind, int *col, int *vecstart, int *vecend, const int *cols)
{
  *vecstart = *col;
  
  while (*col < len && cols[*col] == ind)
    (*col)++;
  
  *vecend = *col - 1;
}



static inline int get_array(int *tmplen, int *current_tmp_size, 
  const int vecstart, const int vecend,
  double *restrict b, int *restrict brows,
  const double *restrict a, const int *restrict rows)
{
  int k;
  
  *tmplen = vecend - vecstart;
  
  if (*tmplen > *current_tmp_size)
  {
    *current_tmp_size = *tmplen;
    
    b = realloc(b, ((*current_tmp_size)+1) * sizeof(*b));
    CHECKMALLOC(b);
    
    brows = realloc(brows, ((*current_tmp_size)+1) * sizeof(*brows));
    CHECKMALLOC(brows);
  }
  
  for (k=0; k<=*tmplen; k++)
  {
    b[k] = a[k + vecstart];
    brows[k] = rows[k + vecstart];
  }
  
  return 0;
}



/**
 * @brief 
 * Compute the cosine similarity matrix of a sparse, COO-stored
 * matrix.
 * 
 * @details
 * The implementation assumes the data is sorted by column index, 
 * i.e. the COO is "column-major".
 * 
 * Note that if the number of rows times the number of columns of
 * the sparse matrix is equal to len, then your matrix is actually
 * dense, but stored in a stupid way.
 * 
 * @param index
 * 0 or 1 indexing from 0 or 1, respectively.
 * @param n
 * The total number of columns of sparsely-stored input matrix x, 
 * i.e., the number of columns of the matrix if it were densely
 * stored.
 * @param len
 * The length of the a/rows/cols vectors.
 * @param a
 * The data for the input matrix, in COO (row, column, value) format.
 * @param rows/cols
 * The row/column index vectors.
 * @param cos
 * The output nxn matrix.
 * 
 * @return
 * The function returns -1 if needed memory cannot be allocated, and
 * 0 otherwise.
*/
int cosine_sparse_coo(const int index, const int n, const int len, 
  const double *restrict a, const int *restrict rows, const int *restrict cols, 
  double *restrict cos)
{
  int i, j, k, l;
  int info;
  int col;
  double xy, xx, yy;
  double tmp;
  
  int vec1start, vec1end;
  int vec2start, vec2end;
  vec1end = 0;
  
  int len_colj;
  int current_tmp_size = TMP_VEC_SIZE;
  double *a_colj = malloc(current_tmp_size * sizeof(*a_colj));
  CHECKMALLOC(a_colj);
  int *rows_colj = malloc(current_tmp_size * sizeof(*rows_colj));
  CHECKMALLOC(rows_colj);
  
  
  set2zero(n*n, cos);
  
  for (j=0; j<n; j++)
  {
    col = vec1end;
    get_startend(len, j+index, &col, &vec1start, &vec1end, cols);
    
    // NaN-out row and column if col is 0
    if (vec1end < vec1start)
    {
      vec1end++;
      set2nan(j, n, cos);
      continue;
    }
    
    // store j't column of data/rows for better cache access
    info = get_array(&len_colj, &current_tmp_size, vec1start, vec1end, a_colj, rows_colj, a, rows);
    if (info) return info;
    
    xx = sparsedot_self(0, len_colj, rows_colj, a_colj);
    xx = 1. / sqrt(xx);
    
    // i'th column, etc.
    for (i=j+1; i<n; i++)
    {
      get_startend(len, i+index, &col, &vec2start, &vec2end, cols);
      
      
      k = 0;
      l = vec2start;
      xy = 0.;
      yy = 0.;
      
      
      while (k <= len_colj && l <= vec2end)
      {
        // catch up row of colj to row of coli
        while (k <= len_colj && rows_colj[k] < rows[l])
          k++;
        
        // dot products
        while (k <= len_colj && l <= vec2end && rows_colj[k] == rows[l])
        {
          tmp = a[l];
          xy += a_colj[k] * tmp;
          yy += tmp*tmp;
          k++;
          l++;
        }
        
        // catch up row of coli to row of colj, self dot product along the way
        if (k <= len_colj)
        {
          while (l <= vec2end && rows_colj[k] > rows[l])
          {
            tmp = a[l];
            yy += tmp*tmp;
            l++;
          }
        }
      }
      
      for (l=l; l<=vec2end; l++)
      {
        tmp = a[l];
        yy += tmp*tmp;
      }
      
      
      if (xy > EPSILON && yy > EPSILON)
        cos[i + n*j] = xy * xx / sqrt(yy);
    }
    
    vec1end++;
  }
  
  
  free(a_colj);
  free(rows_colj);
  
  
  diag2one(n, cos);
  symmetrize(n, cos);
  
  return 0;
}
