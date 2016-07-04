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

#ifndef __COOP_XPOSE_H__
#define __COOP_XPOSE_H__


static inline void xpose(const int m, const int n, const double *const restrict x, double *restrict tx)
{
  const int blocksize = 16; // TODO check cache line explicitly
  
  for (int j=0; j<n; j+=blocksize)
  {
    for (int i=0; i<m; i+=blocksize)
    {
      for (int col=j; col<j+blocksize && col<n; ++col)
      {
        for (int row=i; row<i+blocksize && row<m; ++row)
          tx[col + n*row] = x[row + m*col];
      }
    }
  }
}

/*
void matprinter(int m, int n, double *x)
{
  int i, j;
  
  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
      printf("%.0f ", x[i+m*j]);
    
    putchar('\n');
  }
}

int main()
{
  const int m = 5;
  const int n = 3;
  double *x = malloc(m*n * sizeof(*x));
  double *y = malloc(n*m * sizeof(*y));
  
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      x[i + m*j] = i + m*j + 1;
  
  matprinter(m, n, x);
  putchar('\n');
  
  xpose(m, n, x, y);
  
  matprinter(n, m, y);
  
  // putchar('\n');
  // for (int j=0; j<n; j++)
  //   for (int i=0; i<m; i++)
  //     printf("%f\n", y[i + m*j]);
  
  return 0;
}

*/


#endif
