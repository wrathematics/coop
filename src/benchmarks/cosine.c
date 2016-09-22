#include "../utils/internal/gen.h"
#include "../utils/internal/timer.h"

#include <stdbool.h>

int main()
{
  const int nreps = 100;
  const int m = 10000;
  const int n = 250;
  double *x, *co;
  
  x = gen_runif2(m*n);
  co = zeromat2(n*n);
  
  TIMER(coop_cosine_mat(false, false, m, n, x, co), nreps);
  
  free(x);
  free(co);
  return COOP_OK;
}
