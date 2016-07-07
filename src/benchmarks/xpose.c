#include <stdlib.h>

#include "timer.h"
#include "../utils/safeomp.h"
#include "../utils/xpose.h"

static inline void xpose_naive(const int m, const int n, const double *const restrict x, double *restrict tx)
{
  #pragma omp parallel for default(none) shared(tx) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    for (int i=0; i<m; i++)
      tx[j + n*i] = x[i + mj];
  }
}



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
