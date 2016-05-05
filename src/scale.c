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

#include <stdbool.h>
#include <math.h>
#include "coop.h"
#include "omputils.h"


// center and/or scale x in place
int coop_scale(const bool centerx, const bool scalex, const int m, const int n, double *restrict x, double *restrict colmeans, double *restrict colvars)
{
  int i, j;
  int mj;
  double colmean, colvar;
  double dt, tmp;
  
  if (m == 0 || n == 0) return 0;
  
  // Doing both at once, if needed, is more performant
  if (centerx && scalex)
  {
    tmp = 1. / ((double) m-1);
    
    #pragma omp parallel for private(i, j, colmean, colvar, dt) shared(x) if(m*n > OMP_MIN_SIZE)
    for (j=0; j<n; j++)
    {
      mj = m*j;
      colmean = 0;
      colvar = 0;
      
      SAFE_SIMD
      for (i=0; i<m; i++)
      {
        dt = x[i + mj] - colmean;
        colmean += dt/((double) i+1);
        colvar += dt * (x[i + mj] - colmean);
      }
      
      colvar = sqrt(colvar * tmp);
      
      // Remove mean and variance
      SAFE_SIMD
      for (i=0; i<m; i++)
        x[i + mj] = (x[i + mj] - colmean) / colvar;
      
      colmeans[j] = colmean;
      colvars[j] = colvar;
    }
  }
  else if (centerx)
  {
    const double div = 1. / ((double) m);
    
    #pragma omp parallel for private(i, j, colmean) shared(x) if(m*n > OMP_MIN_SIZE)
    for (j=0; j<n; j++)
    {
      mj = m*j;
      colmean = 0;
      
      // Get column mean
      SAFE_SIMD
      for (i=0; i<m; i++)
        colmean += x[i   + mj] * div;
      
      // Remove mean from column
      SAFE_SIMD
      for (i=0; i<m; i++)
        x[i + mj] -= colmean;
      
      colmeans[j] = colmean;
    }
  }
  else if (scalex) // RMSE
  {
    const double div = 1./((double) m-1);
    
    #pragma omp parallel for private(i, j, colvar, tmp) shared(x) if (m*n > OMP_MIN_SIZE)
    for (j=0; j<n; j++)
    {
      mj = m*j;
      colvar = 0;
      
      // Get column variance
      SAFE_SIMD
      for (i=0; i<m; i++)
      {
        tmp = x[i + mj];
        colvar += tmp*tmp*div;
      }
      
      colvar = sqrt(colvar);
      
      // Remove variance from column
      SAFE_SIMD
      for (i=0; i<m; i++)
        x[i + mj] /= colvar;
      
      colvars[j] = colvar;
    }
  }
  
  return 0;
}
