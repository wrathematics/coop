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


static inline void compute_sums(const int m, const int mi, const double *restrict vec, const double *x, double *restrict sumx, double *restrict sumy, int *restrict len)
{
  int k;
  
  *sumx = 0;
  *sumy = 0;
  *len = 0;
  
  SAFE_SIMD
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



int coop_cosine_mat_inplace_pairwise(const int m, const int n, const double *restrict x, double *restrict cos)
{
  int i, j, k;
  int mj, mi;
  int len;
  double xx, yy, xy, xval, yval;
  double tmp;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (j=0; j<n; j++)
  {
    mj = m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    #pragma omp parallel for private(i, k, mi, xx, yy, xy, xval, yval, len) if(m*n > OMP_MIN_SIZE)
    for (i=j; i<n; i++)
    {
      mi = m*i;
      
      xx = yy = xy = 0.0;
      len = 0;
      
      SAFE_SIMD
      for (k=0; k<m; k++)
      {
        if (!isnan(vec[k]) && !isnan(x[k + mi]))
        {
          xval = vec[k];
          yval = x[k + mi];
          
          xx += xval * xval;
          yy += yval * yval;
          xy += xval * yval;
          len++;
        }
      }
      
      if (len == 0)
      {
        set_na_real(cos + (i + n*j));
        continue;
      }
      
      tmp = xy / sqrt(xx * yy);
      cos[i + n*j] = tmp;
      cos[j + n*i] = tmp;
    }
  }
  
  free(vec);
  
  return 0;
}



int coop_covar_mat_inplace_pairwise(const int m, const int n, const double *restrict x, double *restrict cov)
{
  int i, j, k;
  int mj, mi;
  int len;
  double meanx, meany;
  double mmcp, tmp;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (j=0; j<n; j++)
  {
    mj = m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    #pragma omp parallel for private(i, k, mi, meanx, meany, len, mmcp) if(m*n > OMP_MIN_SIZE)
    for (i=j; i<n; i++)
    {
      mi = m*i;
      
      compute_sums(m, mi, vec, x, &meanx, &meany, &len);
      
      if (len == 0)
      {
        set_na_real(cov + (i + n*j));
        continue;
      }
      
      meanx /= (double) len;
      meany /= (double) len;
      
      mmcp = 0.0;
      SAFE_SIMD
      for (k=0; k<m; k++)
      {
        if (!isnan(vec[k]) && !isnan(x[k + mi]))
          mmcp += (vec[k] - meanx) * (x[k + mi] - meany);
      }
      
      tmp = mmcp * ((double) 1.0/(len-1));
      cov[i + n*j] = tmp;
      cov[j + n*i] = tmp;
    }
  }
  
  free(vec);
  
  return 0;
}
