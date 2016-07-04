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
#include "utils/safeomp.h"
#include "utils/special_vals.h"


static inline void compute_sums(const int m, const int mi, const double * const restrict vec, const double * const restrict x, double *restrict sumx, double *restrict sumy, int *restrict len)
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



int coop_cosine_mat_inplace_pairwise(const int m, const int n, const double * const restrict x, double *restrict cos)
{
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    #pragma omp parallel for default(none) shared(j, vec, cos) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const int mi = m*i;
      
      double xx, xy, yy;
      xx = xy = yy = 0.0;
      int len = 0;
      
      SAFE_SIMD
      for (int k=0; k<m; k++)
      {
        if (!isnan(vec[k]) && !isnan(x[k + mi]))
        {
          const double xval = vec[k];
          const double yval = x[k + mi];
          
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
      
      const double tmp = xy / sqrt(xx * yy);
      cos[i + n*j] = tmp;
      cos[j + n*i] = tmp;
    }
  }
  
  free(vec);
  
  return 0;
}



int coop_pcor_mat_inplace_pairwise(const int m, const int n, const double * const restrict x, double *restrict cor)
{
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    #pragma omp parallel for default(none) shared(j, vec, cor) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const int mi = m*i;
      
      int len;
      double meanx, meany;
      compute_sums(m, mi, vec, x, &meanx, &meany, &len);
      
      if (len == 0 || len == 1)
      {
        set_na_real(cor + (i + n*j));
        set_na_real(cor + (j + n*i));
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
      
      const double tmp = mmcp / sdx / sdy / (dlen - 1.0);
      cor[i + n*j] = tmp;
      cor[j + n*i] = tmp;
    }
  }
  
  free(vec);
  
  return 0;
}



int coop_covar_mat_inplace_pairwise(const int m, const int n, const double * const restrict x, double *restrict cov)
{
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    #pragma omp parallel for default(none) shared(j, vec, cov) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const int mi = m*i;
      
      int len;
      double meanx, meany;
      compute_sums(m, mi, vec, x, &meanx, &meany, &len);
      
      if (len == 0)
      {
        set_na_real(cov + (i + n*j));
        set_na_real(cov + (j + n*i));
        continue;
      }
      
      meanx /= (double) len;
      meany /= (double) len;
      
      double mmcp = 0.0;
      SAFE_SIMD
      for (int k=0; k<m; k++)
      {
        if (!isnan(vec[k]) && !isnan(x[k + mi]))
          mmcp += (vec[k] - meanx) * (x[k + mi] - meany);
      }
      
      const double tmp = mmcp * ((double) 1.0/(len-1));
      cov[i + n*j] = tmp;
      cov[j + n*i] = tmp;
    }
  }
  
  free(vec);
  
  return 0;
}
