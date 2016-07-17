#include <stdlib.h>

#include "../utils/safeomp.h"
#include "../utils/fill.h"

#include "../utils/internal/timer.h"

#include "../utils/naive/symmetrize.h"


int main()
{
  const int nreps = 100;
  const int n = 4000;
  double *x = malloc(n*n * sizeof(*x));
  
  
  TIMER(symmetrize_naive(n, x), nreps);
  TIMER(symmetrize(n, x), nreps);
  
  free(x);
  return 0;
}
