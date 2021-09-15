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


#include <stdlib.h>
#include <string.h>

#include "coop.h"
#include "utils/fill.h"
#include "utils/inverse.h"
#include "utils/safeomp.h"
#include "utils/special_vals.h"

#define TYPE_MIN    1
#define TYPE_MAX    2
#define TYPE_ABSMIN 3
#define TYPE_ABSMAX 4

typedef struct
{
  int K;                // Number of comparisons to retain
  double *restrict co;  // co-variance, sine, ...
  int *restrict I;      // i
  int *restrict J;      // j
  int *restrict L;      // length of run (number of comparisons)
} maxco_t;

// almost sorted sort
// x has some number in the first slot, and is sorted increasing otherwise
// sort x low to high, bring y and z along for the ride
static inline void assort(const register int K, double *restrict A, int *restrict I, int *restrict J, int *restrict L)
{
  int ind;
  const register double atmp = A[0];
  const register int itmp = I[0];
  const register int jtmp = J[0];
  const register int ltmp = L[0];
  
  for (ind=1; ind<K && A[ind] < atmp; ind++)
  {
    A[ind] = A[ind + 1];
    I[ind] = I[ind + 1];
    J[ind] = J[ind + 1];
    L[ind] = L[ind + 1];
  }
  
  A[ind] = atmp;
  I[ind] = itmp;
  J[ind] = jtmp;
  L[ind] = ltmp;
}



// A sorted least to greatest
static inline void rename_me(const int type, const double cmp, const int K, maxco_t *mx)
{
  if (type == TYPE_MIN)
  {
    // TODO
  }
  else if (type == TYPE_ABSMIN)
  {
    // TODO
  }
  else if (type == TYPE_MAX)
  {
    if (cmp > mx->co[0])
      assort(mx->K, mx->co, mx->I, mx->J);
  }
  else if (type == TYPE_ABSMAX)
  {
    // TODO
  }
}




static inline void compute_sums(const int m, const int mi, const double * const restrict vec, const double * const restrict x, double *restrict sumx, double *restrict sumy, int *restrict len)
{
  int k;
  
  *sumx = 0;
  *sumy = 0;
  *len = 0;
  
  PLEASE_VECTORIZE
  for (k=0; k<m; k++)
  {
    if (!isnan(vec[k]) && !isnan(x[k + mi]))
    {
      *sumx += vec[k];
      *sumy += x[k + mi];
      (*len)++;
    }
  }
}



// cor - vals, I/J - their indices, L - length of the run
int coop_maxpcor_mat_inplace_pairwise(const bool inv, const int m, const int n, const double * const restrict x, maxco_t *mx)
{
  int check;
  int ind = 0;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    // #pragma omp parallel for default(none) shared(j, vec) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const int mi = m*i;
      
      int len;
      double meanx, meany;
      compute_sums(m, mi, vec, x, &meanx, &meany, &len);
      
      if (len == 0 || len == 1)
      {
        set_na_real(mx->cor + (i + n*j));
        set_na_real(mx->cor + (j + n*i));
        
        mx->I[ind] = i;
        mx->J[ind] = j;
        mx->L[ind] = len;
        ind++;
        
        continue;
      }
      
      const double dlen = (double) len;
      meanx /= dlen;
      meany /= dlen;
      
      double sdx = 0.;
      double sdy = 0.;
      
      SAFE_SIMD
      for (int k=0; k<m; k++)
      {
        if (!isnan(vec[k]) && !isnan(x[k + mi]))
        {
          sdx += (vec[k] - meanx)*(vec[k] - meanx);
          sdy += (x[k + mi] - meany)*(x[k + mi] - meany);
        }
      }
      
      sdx = sqrt(sdx/(dlen-1.));
      sdy = sqrt(sdy/(dlen-1.));
      
      double mmcp = 0.0;
      SAFE_SIMD
      for (int k=0; k<m; k++)
      {
        if (!isnan(vec[k]) && !isnan(x[k + mi]))
          mmcp += (vec[k] - meanx) * (x[k + mi] - meany);
      }
      
      rename_me(TYPE_MAX, mmcp*dlen, mx);
      // cor[i + n*j] = mmcp / sdx / sdy / (dlen - 1.0);;
    }
  }
  
  free(vec);
  
  if (inv)
  {
    check = inv_sym_chol(n, cor);
    CHECKRET(check);
  }
  
  
  return COOP_OK;
}
