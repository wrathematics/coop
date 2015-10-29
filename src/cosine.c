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
#include <math.h>
#include "cosine.h"
#include "omp.h"


void dgemm_(const char *transa, const char *transb, const int *m, const int *n, 
            const int *k, const double *restrict alpha, const double *restrict a, 
            const int *lda, const double *restrict b, const int *ldb, 
            const double *beta, double *restrict c, const int *ldc);


static inline double ddot(const int n, const double *restrict x, const double *restrict y)
{
  int one = 1;
  double dot;
  
  dgemm_(&(char){'t'}, &(char){'n'}, &one, &one, &n, 
    &(double){1.0}, x, &n, y, &n, &(double){0.0}, &dot, &one);
  
  return dot;
}


void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k, 
            const double *restrict alpha, const double *restrict a, const int *lda, 
            const double *restrict beta, double *restrict c, const int *ldc);



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



// Copy upper triangle to lower
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




#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


static double sparsedot(const int vec1start, const int vec1end, 
                        const int vec2start, const int vec2end,
                        const int *rows, const double *a)
{
  const int len = MIN(vec1end-vec1start, vec2end-vec2start);
  int col1 = 0;
  int col2 = 0;
  
  double dot = 0.0;
  
  while (col1 < len && col2 < len)
  {
    while (rows[col1] < rows[col2])
      col1++;
    
    while (rows[col1] == rows[col2])
    {
      dot += a[col1] * a[col2];
      col1++;
      col2++;
    }
    
    while (rows[col1] > rows[col2])
      col2++;
  }
}



static inline double sparsedot_self(const int vecstart, const int vecend, const int *rows, const double *a)
{
  int i = vecstart;
  double dot = 0.0;
  
  while (i <= vecend)
    dot += a[i]*a[i];
  
  return dot;
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
*/
void cosine_sparse_coo(const int n, const int len, const double *restrict a, const int *restrict rows, const int *restrict cols, double *restrict cos)
{
  // TODO note, assuming sorted by column index, then row
  int i, j;
  int row, col;
  double xy, xx, yy;
  
  int vec1start, vec1end;
  int vec2start, vec2end;
  
  col = 0;
  
  for (j=0; j<len; j++)
  {
    vec1start = col;
    while (cols[col] == j) //FIXME 0/1 indexing
      col++;
    vec1end = col - 1;
    
    if (vec1end - vec1start == 0)
      continue;
    
    for (i=0; i<len; i++)
    {
      vec2start = col;
      while (cols[col] == i)
        col++;
      vec2end = col - 1;
      
      xy = sparsedot(vec1start, vec1end, vec2start, vec2end, rows, a);
      if (xy > 1e-8)
      {
        xx = sparsedot_self(vec1start, vec1end, rows, a);
        yy = sparsedot_self(vec2start, vec2end, rows, a);
        
        xy /= sqrt(xx * yy);
      }
      
      cos[i + n*j] = xy;
    }
  }
  
  diag2one(n, cos);
  symmetrize(n, cos);
}


