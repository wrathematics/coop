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

// Functions for computing cosine similarity on sparse data inputs
// TODO: covariance and (pearson) correlation

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "coop.h"
#include "omputils.h"

#define TMP_VEC_SIZE 1024


// ---------------------------------------------
//  Static utils
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
  void *realloc_ptr;

  *tmplen = vecend - vecstart;

  if (*tmplen > *current_tmp_size)
  {
    *current_tmp_size = *tmplen;

    realloc_ptr = realloc(b, ((*current_tmp_size)+1) * sizeof(*b));
    if (realloc_ptr == NULL)
    {
      free(b);
      free(brows);
      return -1;
    }
    else
      b = realloc_ptr;

    realloc_ptr = realloc(brows, ((*current_tmp_size)+1) * sizeof(*brows));
    if (realloc_ptr == NULL)
    {
      free(b);
      free(brows);
      return -1;
    }
    else
      brows = realloc_ptr;
  }

  for (k=0; k<=*tmplen; k++)
  {
    b[k] = a[k + vecstart];
    brows[k] = rows[k + vecstart];
  }

  return 0;
}





// ---------------------------------------------
//  Cosine
// ---------------------------------------------

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
int coop_cosine_sparse_coo(const int index, const int n, const int len,
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

  coop_diag2one(n, cos);
  coop_symmetrize(n, cos);

  return 0;
}
