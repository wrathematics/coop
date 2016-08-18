#include "../utils/safeomp.h"
#include "../utils/fill.h"

#include "../utils/internal/gen.h"
#include "../utils/internal/timer.h"

#include "../utils/naive/symmetrize.h"


int main()
{
  const int nreps = 100;
  const int n = 5000;
  double *x;
  
  x = gen_runif2(n*n);
  
  TIMER(symmetrize_naive(n, x), nreps);
  TIMER(symmetrize(n, x), nreps);
  
  return 0;
}
