#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../utils/safeomp.h"
#include "../utils/xpose.h"

#include "../tests/naive.h"


int main()
{
  const int m = 10;
  const int n = 3;
  double *x = malloc(m*n * sizeof(*x));
  double *tx = malloc(n*m * sizeof(*tx));
  double *truth = malloc(n*m * sizeof(*truth));
  
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      x[i + m*j] = i + m*j + 1;
  
  memset(truth, 0.0, m*n*sizeof(*truth));
  xpose_naive(m, n, x, truth);
  
  memset(tx, 0.0, m*n*sizeof(*tx));
  xpose(m, n, x, tx);
  
  for (int i=0; i<m*n; i++)
  {
    if (truth[i] != tx[i])
    {
      printf("FAIL\n");
      
      matprinter(n, m, tx);
      matprinter(n, m, truth);
      return -1;
    }
  }
  
  printf("PASS\n");
  
  free(x);
  free(tx);
  return 0;
}
