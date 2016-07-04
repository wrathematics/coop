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


#include <math.h>
#include "coop.h"
#include "utils/safeomp.h"


// set diagonal of nxn matrix x to 1
void coop_diag2one(const unsigned int n, double *restrict x)
{
  int i;
  
  SAFE_FOR_SIMD
  for (i=0; i<n; i++)
    x[i + n*i] = 1.0;
}



// Copy lower triangle to upper
void coop_symmetrize(const int n, double *restrict x)
{
  #pragma omp parallel for default(none) shared(x) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    const int nj = n*j;
    
    SAFE_SIMD
    for (int i=j+1; i<n; i++)
      x[j + n*i] = x[i + nj];
  }
}



// replaces upper triangle of the crossproduct of a matrix with its cosine similarity
void coop_fill(const unsigned int n, double *restrict cp)
{
  #pragma omp parallel for default(none) shared(cp) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    const double diagj = cp[j + n*j];
    
    const int nj = n*j;
    cp[j + nj] = 1;
    
    SAFE_SIMD
    for (int i=j+1; i<n; i++)
      cp[i + nj] /= sqrt(cp[i + n*i] * diagj);
  }
}
