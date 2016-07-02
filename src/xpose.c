#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "omputils.h"


void xpose(const int m, const int n, const double *const restrict x, double *restrict tx)
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
