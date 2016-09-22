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


#ifndef __COOP_LIB_SCALE_H__
#define __COOP_LIB_SCALE_H__


static inline void centerscale(const int j, const int m, const int n, double *restrict x, double *restrict colmean, double *restrict colvar)
{
  const double tmp = 1. / ((double) m-1);
  
  const int mj = m*j;
  *colmean = 0;
  *colvar = 0;
  
  SAFE_FOR_SIMD
  for (int i=0; i<m; i++)
  {
    double dt = x[i + mj] - *colmean;
    *colmean += dt/((double) i+1);
    *colvar += dt * (x[i + mj] - *colmean);
  }
  
  *colvar = sqrt(*colvar * tmp);
  
  // Remove mean and variance
  SAFE_FOR_SIMD
  for (int i=0; i<m; i++)
    x[i + mj] = (x[i + mj] - *colmean) / *colvar;
}



static inline double center(const int j, const int m, const int n, double *x)
{
  const double div = 1. / ((double) m);
  
  const int mj = m*j;
  double colmean = 0;
  
  // Get column mean
  SAFE_FOR_SIMD
  for (int i=0; i<m; i++)
    colmean += x[i + mj] * div;
  
  // Remove mean from column
  SAFE_FOR_SIMD
  for (int i=0; i<m; i++)
    x[i + mj] -= colmean;
  
  return colmean;
}



static inline double scale(const int j, const int m, const int n, double *x)
{
  const double div = 1./((double) m-1);
  
  const int mj = m*j;
  double colvar = 0;
  
  // Get column variance
  SAFE_FOR_SIMD
  for (int i=0; i<m; i++)
  {
    double tmp = x[i + mj];
    colvar += tmp*tmp*div;
  }
  
  colvar = sqrt(colvar);
  
  // Remove variance from column
  SAFE_FOR_SIMD
  for (int i=0; i<m; i++)
    x[i + mj] /= colvar;
  
  return colvar;
}



static inline int scale_nostore(const bool centerx, const bool scalex, const int m, const int n, double *restrict x)
{
  if (m == 0 || n == 0)
    return COOP_OK;
  
  // Doing both at once, if needed, is more performant
  if (centerx && scalex)
  {
    double colmean;
    double colvar;
    #pragma omp parallel for shared(x) if (m*n > OMP_MIN_SIZE)
    for (int j=0; j<n; j++)
      centerscale(j, m, n, x, &colmean, &colvar);
    
  }
  else if (centerx)
  {
    #pragma omp parallel for shared(x) if (m*n > OMP_MIN_SIZE)
    for (int j=0; j<n; j++)
      center(j, m, n, x);
  }
  else if (scalex) // RMSE
  {
    #pragma omp parallel for shared(x) if (m*n > OMP_MIN_SIZE)
    for (int j=0; j<n; j++)
      scale(j, m, n, x);
  }
  
  return COOP_OK;
}


#endif
