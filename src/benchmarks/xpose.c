#include <stdlib.h>

#include "../utils/safeomp.h"
#include "../utils/xpose.h"

#include "../utils/internal/timer.h"

#include "../utils/naive/xpose.h"


int main()
{
  const int nreps = 100;
  const int m = 10000;
  const int n = 250;
  double *x = malloc(m*n * sizeof(*x));
  double *tx = malloc(n*m * sizeof(*tx));
  
  TIMER(xpose_naive(m, n, x, tx), nreps);
  TIMER(xpose(m, n, x, tx), nreps);
  
  free(x);
  free(tx);
  return 0;
}
