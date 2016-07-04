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

#ifndef __COOP_SUMSTATS_H__
#define __COOP_SUMSTATS_H__


// x[*, j] -= colmean(x[*, j])
static void remove_colmeans(const int m, const int n, double *restrict x)
{
  if (m == 0 || n == 0)
    return;
    
  const double div = 1. / ((double) m);
  
  #pragma omp parallel for default(none) shared(x) if(m*n > OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    double colmean = 0;
    
    // Get column mean
    SAFE_SIMD
    for (int i=0; i<m; i++)
      colmean += x[i + m*j];
      
    colmean *= div;
    
    // Remove mean from column
    SAFE_SIMD
    for (int i=0; i<m; i++)
      x[i + m*j] -= colmean;
  }
}



// same as above but return the means vector
static void remove_colmeans_retmean(const int m, const int n, double *restrict x, double *restrict colmeans)
{
  if (m == 0 || n == 0)
    return;
    
  const double div = 1. / ((double) m);
  
  #pragma omp parallel for default(none) shared(x, colmeans) if(m*n > OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    colmeans[j] = 0;
    
    // Get column mean
    SAFE_SIMD
    for (int i=0; i<m; i++)
      colmeans[j] += x[i + m*j];
      
    colmeans[j] *= div;
    
    // Remove mean from column
    SAFE_SIMD
    for (int i=0; i<m; i++)
      x[i + m*j] -= colmeans[j];
  }
}



// compute the mean of a vector
static inline double mean(const int n, const double * const restrict x)
{
  const double divbyn = 1. / ((double) n);
  double mean = 0.;
  
  SAFE_FOR_SIMD
  for (int i=0; i<n; i++)
    mean += x[i];
  
  return mean*divbyn;
}



#endif
