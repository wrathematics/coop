#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../utils/safeomp.h"
#include "../utils/xpose.h"

#include "../utils/naive/all_equal.h"
#include "../utils/naive/xpose.h"


int main()
{
  int ret;
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
  
  ret = all_equal("Transpose", true, n, m, tx, truth);
  
  free(x);
  free(tx);
  free(truth);
  return ret;
}
