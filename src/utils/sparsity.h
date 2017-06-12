/*  Copyright (c) 2015-2016 Drew Schmidt
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

#ifndef __COOP_LIB_SPARSITY_H__
#define __COOP_LIB_SPARSITY_H__


// Number of 0's for integer matrix
static int coop_sparsity_int(const int m, const int n, const int * const restrict x)
{
  int count = 0;
  
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    
    for (int i=0; i<m; i++)
    {
      if (x[i + mj] == 0)
        count++;
    }
  }
  
  return count;
}



// Number of (approximate) 0's for double matrix
static int coop_sparsity_dbl(const int m , const int n, const double * const restrict x, const double tol)
{
  int count = 0;

  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    
    for (int i=0; i<m; i++)
    {
      if (fabs(x[i + mj]) < tol)
        count++;
    }
  }
  
  return count;
}


#endif
