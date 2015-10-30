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
  int one = 1;
  double dot;
  
  dgemm_(&(char){'t'}, &(char){'n'}, &one, &one, &n, 
    &(double){1.0}, x, &n, y, &n, &(double){0.0}, &dot, &one);
  
  return dot;
}



// upper triangle of t(x) %*% x
static inline int crossprod(const int m, const int n, const double *restrict x, const double alpha, double *restrict c)
{
  int info = 0;
  
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &alpha, x, &m, &(double){0.0}, c, &n);
  
  return info;
}



static inline void diag2one(const unsigned int n, double *restrict x)
{
  int i;
  
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
    x[i + n*i] = 1.0;
}



// replaces upper triangle of the crossproduct of a matrix with its cosine similarity
static inline void fill(const unsigned int n, double *restrict crossprod)
{
  int i, j;
  double diagj;
  
  #pragma omp parallel for private(i,j,diagj) default(shared) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    diagj = crossprod[j + n*j];
    
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      crossprod[i + n*j] /= sqrt(crossprod[i + n*i] * diagj);
  }
  
  diag2one(n, crossprod);
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
  crossprod(m, n, x, 1.0, cos);
  fill(n, cos);
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
  
  double cp = ddot(n, x, y);
  
  crossprod(n, 1, x, 1.0, &normx);
  crossprod(n, 1, y, 1.0, &normy);
  
  return cp / sqrt(normx * normy);
}





// ---------------------------------------------
//  Sparse
// ---------------------------------------------

static inline void set2zero(const unsigned int n, double *restrict x)
{
  memset(x, 0.0, n*sizeof(*x));
}



static double sparsedot(const int vec1start, const int vec1end, 
                        const int *restrict rows1, const double *restrict a1,
                        const int vec2start, const int vec2end,
                        const int *restrict rows2, const double *restrict a2)
{
  int vec1 = vec1start;
  int vec2 = vec2start;
  
  double dot = 0.0;
  
  
  while (vec1 <= vec1end && vec2 <= vec2end)
  {
    while (rows1[vec1] < rows2[vec2] && vec1 <= vec1end)
      vec1++;
    
    while (rows1[vec1] == rows2[vec2] && vec1 <= vec1end && vec2 <= vec2end)
    {
      dot += a1[vec1] * a2[vec2];
      vec1++;
      vec2++;
    }
    
    while (rows1[vec1] > rows2[vec2] && vec2 <= vec2end)
      vec2++;
  }
  
  return dot;
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
static inline void get_startend(int i, int *col, int *vecstart, int *vecend, const int *cols)
{
  // FIXME 0/1 indexing
  *vecstart = *col;
  
  while (cols[*col] == i)
    (*col)++;
  
  *vecend = *col - 1;
}



/**
 * @brief 
 * Compute the cosine similarity matrix of a sparse, COO-stored
 * matrix.
 * 
 * @details
 * TODO
 * 
 * Note that if the number of rows times the number of columns of
 * the sparse matrix is equal to len, then your matrix is actually
 * dense, but stored in a stupid way.
 * 
 * @param n
 * The total number of columns of sparsely-stored input matrix x, 
 * i.e., the number of columns of the matrix if it were densely
 * stored.
 * @param len
 * The length of the a/rows/cols vectors.
 * @param a
 * The input matrix, in COO (row, column, value) format.
 * @param rows/cols
 * The row/column index vectors.
 * @param cos
 * The output nxn matrix.
 * 
 * @return
 * The function returns -1 if needed memory cannot be allocated, and
 * 0 otherwise.
*/
int cosine_sparse_coo(const int n, const int len, const double *restrict a, const int *restrict rows, const int *restrict cols, double *restrict cos)
{
  // TODO note, assuming sorted by column index, then row
  int i, j, k;
  int row, col;
  double xy, xx, yy;
  double tmp;
  
  int vec1start, vec1end;
  int vec2start, vec2end;
  vec1end = 0;
  
  int tmpend;
  int current_tmp_size = TMP_VEC_SIZE;
  double *tmpa = malloc(current_tmp_size * sizeof(*tmpa));
  checkmalloc(tmpa);
  int *tmprows = malloc(current_tmp_size * sizeof(*tmprows));
  checkmalloc(tmprows);
  
  
  set2zero(n*n, cos);
  
  for (j=0; j<n; j++)
  {
    col = vec1end;
    get_startend(j, &col, &vec1start, &vec1end, cols);
    
    // NaN-out row and column if col is 0
    if (vec1end < vec1start)
    {
      vec1end++;
      
      for (i=j; i<n; i++)
        cos[i + n*j] = NAN;
      
      for (i=0; i<j; i++)
        cos[j + n*i] = NAN;
      
      continue;
    }
    
    // store j't column of data/rows for better cache access
    tmpend = vec1end - vec1start;
    
    if (tmpend > current_tmp_size)
    {
      //TODO check
      tmpa = realloc(tmpa, tmpend * sizeof(*tmpa));
      checkmalloc(tmpa);
      tmprows = realloc(tmprows, tmpend * sizeof(*tmprows));
      checkmalloc(tmprows);
    }
    
    for (k=0; k<=tmpend; k++)
    {
      tmpa[k] = a[k + vec1start];
      tmprows[k] = rows[k + vec1start];
    }
    
    
    // i'th column, etc.
    for (i=j+1; i<n; i++)
    {
      get_startend(i, &col, &vec2start, &vec2end, cols);
      
      xy = sparsedot(0, tmpend, tmprows, tmpa, vec2start, vec2end, rows, a);
      
      if (xy > 1e-10)
      {
        xx = sparsedot_self(0, tmpend, tmprows, tmpa);
        yy = sparsedot_self(vec2start, vec2end, rows, a);
        
        tmp = sqrt(xx * yy);
        if (tmp > 0)
          cos[i + n*j] = xy/tmp;
      }
    }
    
    vec1end++;
  }
  
  
  free(tmpa);
  free(tmprows);
  
  
  diag2one(n, cos);
  symmetrize(n, cos);
  
  return 0;
}


