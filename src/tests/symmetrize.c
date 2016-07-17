#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../utils/safeomp.h"
#include "../utils/fill.h"

#include "../utils/naive/all_equal.h"
#include "../utils/naive/symmetrize.h"

int main()
{
  int ret;
  const int n = 10;
  double *x = malloc(n*n * sizeof(*x));
  double *sym = malloc(n*n * sizeof(*sym));
  double *truth = malloc(n*n * sizeof(*truth));
  
  for (int j=0; j<n; j++)
    for (int i=0; i<n; i++)
      x[i + n*j] = i + n*j + 1;
  
  memcpy(truth, x, n*n*sizeof(*truth));
  symmetrize_naive(n, truth);
  
  memcpy(sym, x, n*n*sizeof(*sym));
  symmetrize(n, sym);
  
  ret = all_equal("Symmetrize", true, n, n, sym, truth);
  
  free(x);
  free(sym);
  free(truth);
  return 0;
}
