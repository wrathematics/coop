#include <stdlib.h>
#include <stdio.h>

#include "timer.h"
#include "../utils/safeomp.h"
#include "../utils/fill.h"

#include "../tests/naive.h"


int main()
{
  const int nreps = 100;
  const int n = 2000;
  double *x = malloc(n*n * sizeof(*x));
  
  
  TIMER(symmetrize_naive(n, x), nreps);
  TIMER(symmetrize(n, x), nreps);
  
  free(x);
  return 0;
}
