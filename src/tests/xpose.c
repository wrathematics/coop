#include "../utils/safeomp.h"
#include "../utils/xpose.h"

#include "../utils/internal/all_equal.h"
#include "../utils/internal/gen.h"

#include "../utils/naive/xpose.h"


int main()
{
  int ret;
  const int m = 10;
  const int n = 3;
  double *x, *tx, *truth;
  
  printf("# Transpose: ");
  
  x = gen_runif2(m*n);
  
  truth = zeromat2(m*n);
  xpose_naive(m, n, x, truth);
  
  tx = zeromat2(m*n);
  xpose(m, n, x, tx);
  
  ret = all_equal(true, n, m, tx, truth);
  
  
  free(x);
  free(tx);
  free(truth);
  return ret;
}
