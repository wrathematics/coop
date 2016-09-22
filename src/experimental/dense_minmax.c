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


#include <stdlib.h>
#include <string.h>

#include "coop.h"
#include "utils/safeomp.h"

#define TYPE_MIN    1
#define TYPE_MAX    2
#define TYPE_ABSMIN 3
#define TYPE_ABSMAX 4


// almost sorted sort
// x has some number in the first slot, and is sorted increasing otherwise
// sort x low to high, bring y and z along for the ride
static inline void assort(const register int k, register double *restrict A, register int *restrict I, register int *restrict J)
{
  int ind;
  const register double atmp = A[0];
  const register int itmp = I[0];
  const register int jtmp = J[0];
  
  for (ind=1; ind<k && A[ind] < atmp; ind++)
  {
    A[ind] = A[ind + 1];
    I[ind] = I[ind + 1];
    J[ind] = J[ind + 1];
  }
  
  A[ind] = atmp;
  I[ind] = itmp;
  J[ind] = jtmp;
}



// A sorted least to greatest
static inline void rename_me(const int type, const double mmcp, const double denom_cov, const int k, double *restrict A, int *restrict I, int *restrict J)
{
  const double tmp = mmcp * denom_cov;
  
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
    if (tmp > A[0])
      assort(k, A, I, J);
  }
  else if (type == TYPE_ABSMAX)
  {
    // TODO
  }
}



// O(m+n) storage
static int co_mat_minmax(const int type, const int m, const int n, const double * const restrict x, const int k, double *restrict A, int *restrict I, int *restrict J)
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
  #pragma omp parallel for default(none) shared(means) if (m*n > OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    
    means[j] = 0.0;
    SAFE_SIMD
    for (int i=0; i<m; i++)
      means[j] += x[i + mj];
    
    means[j] *= denom_mean;
  }
  
  
  // co-operation
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    const double meanx = means[j];
    PLEASE_VECTORIZE
    for (int k=0; k<m; k++)
      vec[k] -= meanx;
    
    #pragma omp parallel for default(none) shared(j, means, vec, A, I, J) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const int mi = m*i;
      
      const double meany = means[i];
      
      double mmcp = 0.0;
      SAFE_SIMD
      for (int l=0; l<m; l++)
        mmcp += vec[l] * (x[l + mi] - meany);
      
      
      // pick top K
      rename_me(type, mmcp, denom_cov, k, A, I, J);
    }
  }
  
  free(vec);
  free(means);
  
  return COOP_OK;
}
