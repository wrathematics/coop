/*  Copyright (c) 2016 Drew Schmidt
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

#include "utils/safeomp.h"

#include "coop.h"
#include "utils/fill.h"
#include "utils/inverse.h"



/*
// O(1) storage
static int coop_covar_vecvec_inplace(const int n, const double * const restrict x, const double * const restrict y, double *restrict cor)
{
  int i;
  const double denom = (double) 1/n;
  double meanx;
  double meany;       // :DDD
  double mmcp = 0.0;  // minus-mean-crossproduct
  
  meanx = 0.0;
  meany = 0.0;
  
  PLEASE_VECTORIZE
  for (i=0; i<n; i++)
  {
    meanx += x[i];
    meany += y[i];
  }
  
  meanx *= denom;
  meany *= denom;
  
  PLEASE_VECTORIZE
  for (i=0; i<n; i++)
    mmcp += (x[i] - meanx) * (y[i] - meany);
    
  *cor = mmcp / ((double)(n-1));
  
  return COOP_OK;
}
*/



// O(m+n) storage
static int co_mat_inplace(const int m, const int n, const double * const restrict x, double *restrict cov)
{
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  double *means = malloc(n * sizeof(*means));
  if (means==NULL)
  {
    free(vec);
    return -1;
  }
  const double denom_mean = (double) 1./m;
  const double denom_cov = (double) 1./(m-1);
  
  
  // get column means
  #pragma omp parallel for shared(means) if (m*n > OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    const size_t mj = m*j;
    
    means[j] = 0.0;
    SAFE_SIMD
    for (int i=0; i<m; i++)
      means[j] += x[i + mj];
    
    means[j] *= denom_mean;
  }
  
  
  // co-operation
  for (int j=0; j<n; j++)
  {
    const size_t mj = m*j;
    
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    const double meanx = means[j];
    PLEASE_VECTORIZE
    for (int k=0; k<m; k++)
      vec[k] -= meanx;
    
    #pragma omp parallel for shared(j, means, vec, cov) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const size_t mi = m*i;
      const double meany = means[i];
      
      double mmcp = 0.0;
      SAFE_SIMD
      for (int k=0; k<m; k++)
        mmcp += vec[k] * (x[k + mi] - meany);
        
      cov[i + n*j] = mmcp * denom_cov;
    }
  }
  
  free(vec);
  free(means);
  
  return COOP_OK;
}



// ---------------------------------------------
//  Interface
// ---------------------------------------------

int coop_pcor_mat_inplace(const bool inv, const int m, const int n, const double * const restrict x, double *restrict cor)
{
  int check = co_mat_inplace(m, n, x, cor);
  CHECKRET(check);
  
  cosim_fill(n, cor);
  
  if (inv)
  {
    check = inv_sym_chol(n, cor);
    CHECKRET(check);
  }
  
  symmetrize(n, cor);
  
  return COOP_OK;
}



int coop_covar_mat_inplace(const bool inv, const int m, const int n, const double * const restrict x, double *restrict cov)
{
  int check = co_mat_inplace(m, n, x, cov);
  CHECKRET(check);
  
  if (inv)
  {
    check = inv_sym_chol(n, cov);
    CHECKRET(check);
  }
  
  symmetrize(n, cov);
  
  return COOP_OK;
}
