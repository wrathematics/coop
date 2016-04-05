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

// Functions for computing covariance, (pearson) correlation, and cosine similarity

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "coop.h"
#include "omp.h"


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
int coop_covar_vecvec_inplace(const int n, const double *restrict x, const double *restrict y, double *restrict cor)
{
  int i;
  const double denom = (double) 1/n;
  double meanx;
  double meany;       // :DDD
  double mmcp = 0.0;  // minus-mean-crossproduct
  
  meanx = 0.0;
  meany = 0.0;
  
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
  {
    meanx += x[i];
    meany += y[i];
  }
  
  meanx *= denom;
  meany *= denom;
  
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
    mmcp += (x[i] - meanx) * (y[i] - meany);
  
  *cor = mmcp / ((double)(n-1));
  
  return 0;
}



// storage: m+n doubles
int coop_covar_mat_inplace(const int m, const int n, const double *restrict x, double *restrict cov)
{
  int i, j, k;
  int mj, mi;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  double *means = malloc(n * sizeof(*means));
  if (x==NULL)
  {
    
  }
  
  const double denom_mean = (double) 1./m;
  const double denom_cov = (double) 1./(m-1);
  double meanx;
  double meany; // :DDD
  double mmcp;  // minus-mean-crossproduct
  
  for (j=0; j<n; j++)
  {
    mj = m*j;
    
    means[j] = 0.0;
    SAFE_FOR_SIMD
    for (i=0; i<m; i++)
      means[j] += x[i + mj];
    
    means[j] *= denom_mean;
  }
  
  
  for (j=0; j<n; j++)
  {
    mj = m*j;
    
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    meanx = means[j];
    SAFE_FOR_SIMD
    for (k=0; k<m; k++)
      vec[k] -= meanx;
    
    #pragma omp parallel for private(i, mi, meany, mmcp)
    for (i=j; i<n; i++)
    {
      mi = m*i;
      
      meany = means[i];
      
      mmcp = 0.0;
      SAFE_SIMD
      for (k=0; k<m; k++)
        mmcp += vec[k] * (x[k + mi] - meany);
      
      cov[i + n*j] = mmcp * denom_cov;
    }
  }
  
  coop_symmetrize(n, cov);
  free(vec);
  free(means);
  
  return 0;
}
