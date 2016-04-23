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
#include "omp.h"


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
  int i, j;
  int nj;
  
  #pragma omp parallel for private(i) default(shared) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    nj = n*j;
    
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      x[j + n*i] = x[i + nj];
  }
}



// replaces upper triangle of the crossproduct of a matrix with its cosine similarity
void cosim_fill(const unsigned int n, double *restrict cp)
{
  int i, j;
  int nj;
  double diagj;
  
  #pragma omp parallel for private(i,j,diagj) default(shared) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    diagj = cp[j + n*j];
    
    nj = n*j;
    cp[j + nj] = 1;
    
    SAFE_SIMD
    for (i=j+1; i<n; i++)
      cp[i + nj] /= sqrt(cp[i + n*i] * diagj);
  }
}



// Number of 0's for integer matrix
int coop_sparsity_int(const int m, const int n, const int *x)
{
  int i, j;
  int mj;
  int count = 0;
  
  for (j=0; j<n; j++)
  {
    mj = m*j;
    
    for (i=0; i<m; i++)
    {
      if (x[i + mj] == 0)
        count++;
    }
  }
  
  return count;
}



// Number of (approximate) 0's for double matrix
int coop_sparsity_dbl(const int m , const int n, double *x, const double tol)
{
  int i, j;
  int mj;
  int count = 0;
  
  for (j=0; j<n; j++)
  {
    mj = m*j;
    
    for (i=0; i<m; i++)
    {
      if (fabs(x[i + mj]) < tol)
        count++;
    }
  }
  
  return count;
}
