#include "../utils/safeomp.h"
#include "../utils/fill.h"

#include "../utils/internal/all_equal.h"
#include "../utils/internal/gen.h"

#include "../utils/naive/symmetrize.h"

int main()
{
  int ret;
  const int n = 10;
  double *x, *sym, *truth;
  
  printf("# Symmetrize: ");
  
  x = gen_boring2(n*n);
  
  truth = cpalloc2(n*n, x);
  symmetrize_naive(n, truth);
  
  sym = cpalloc2(n*n, x);
  symmetrize(n, sym);
  
  ret = all_equal(true, n, n, sym, truth);
  
  
  free(x);
  free(sym);
  free(truth);
  return ret;
}
