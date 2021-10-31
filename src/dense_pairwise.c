/*  Copyright (c) 2016, 2021 Drew Schmidt
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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils/safeomp.h"

#include "coop.h"
#include "utils/fill.h"
#include "utils/inverse.h"
#include "utils/special_vals.h"


static inline void compute_sums(const int m, const size_t mi,
  const double *const restrict vec, const double *const restrict x,
  double *restrict sumx, double *restrict sumy, int *restrict len)
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



// -----------------------------------------------------------------------------
// cosine
// -----------------------------------------------------------------------------
int coop_cosine_mat_inplace_pairwise(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cos)
{
  int check;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<n; j++)
  {
    const size_t mj = (size_t)m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    const size_t nj = (size_t)n*j;
    
    #pragma omp parallel for shared(j, vec, cos) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const size_t mi = (size_t)m*i;
      
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
        set_na_real(cos + (i + nj));
        continue;
      }
      
      cos[i + nj] = xy / sqrt(xx * yy);
    }
  }
  
  free(vec);
  
  if (inv)
  {
    check = inv_sym_chol(n, cos);
    CHECKRET(check);
  }
  
  symmetrize(n, cos);
  
  return COOP_OK;
}



// -----------------------------------------------------------------------------
// correlation
// -----------------------------------------------------------------------------
int coop_pcor_mat_inplace_pairwise(const bool inv, const int m, const int n, const double * const restrict x, double *restrict cor)
{
  int check;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<n; j++)
  {
    const size_t mj = (size_t)m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    const size_t nj = (size_t)n*j;
    
    #pragma omp parallel for shared(j, vec, cor) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const size_t mi = (size_t)m*i;
      
      int len;
      double meanx, meany;
      compute_sums(m, mi, vec, x, &meanx, &meany, &len);
      
      if (len == 0 || len == 1)
      {
        set_na_real(cor + (i + nj));
        set_na_real(cor + (j + (size_t)n*i));
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
      
      cor[i + nj] = mmcp / sdx / sdy / (dlen - 1.0);;
    }
  }
  
  free(vec);
  
  if (inv)
  {
    check = inv_sym_chol(n, cor);
    CHECKRET(check);
  }
  
  symmetrize(n, cor);
  
  return COOP_OK;
}



// -----------------------------------------------------------------------------
// covariance
// -----------------------------------------------------------------------------
int coop_covar_mat_inplace_pairwise(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cov)
{
  int check;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<n; j++)
  {
    const size_t mj = (size_t)m*j;
    memcpy(vec, x+mj, m*sizeof(*vec));
    
    const size_t nj = (size_t)n*j;
    
    #pragma omp parallel for shared(j, vec, cov) if(m*n > OMP_MIN_SIZE)
    for (int i=j; i<n; i++)
    {
      const size_t mi = (size_t)m*i;
      
      int len;
      double meanx, meany;
      compute_sums(m, mi, vec, x, &meanx, &meany, &len);
      
      if (len == 0)
      {
        set_na_real(cov + (i + nj));
        set_na_real(cov + (j + (size_t)n*i));
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
      
      cov[i + nj] = mmcp * ((double) 1.0/(len-1));
    }
  }
  
  free(vec);
  
  if (inv)
  {
    check = inv_sym_chol(n, cov);
    CHECKRET(check);
  }
  
  symmetrize(n, cov);
  
  return COOP_OK;
}



int coop_covar_matmat_inplace_pairwise(const bool inv, const int m, const int nx,
  const double *const restrict x, const int ny, const double *const restrict y,
  double *restrict cov)
{
  int check;
  double *vec = malloc(m * sizeof(*vec));
  CHECKMALLOC(vec);
  
  
  for (int j=0; j<ny; j++)
  {
    const size_t mj = (size_t)m*j;
    memcpy(vec, y+mj, m*sizeof(*vec));
    
    
    #pragma omp parallel for shared(j, vec, cov) if(m*nx > OMP_MIN_SIZE)
    for (int i=0; i<nx; i++)
    {
      const size_t mi = (size_t)m*i;
      
      int len;
      double meanx, meany;
      compute_sums(m, mi, vec, x, &meanx, &meany, &len);
      
      if (len == 0)
      {
        set_na_real(cov + (i + nx*j));
        set_na_real(cov + (j + nx*i));
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
      
      cov[i + nx*j] = mmcp * ((double) 1.0/(len-1));
    }
  }
  
  free(vec);
  
  if (nx == ny && inv)
  {
    check = inv_gen_lu(nx, cov);
    CHECKRET(check);
  }
  
  return COOP_OK;
}
