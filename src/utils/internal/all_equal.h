#ifndef __COOP_UTILS_INTERNAL_ALLEQUAL_H__
#define __COOP_UTILS_INTERNAL_ALLEQUAL_H__

#include <stdbool.h>
#include <math.h>
#include <stdio.h>

#include "printer.h"

#define EPS 1e-10
#define CMP(a,b) (fabs(a-b) < EPS)

static inline int all_equal(const char *testname, const bool printonfail, const int m, const int n, const double *restrict const test, const double *restrict const truth)
{
  printf("### %s: ", testname);
  
  for (int i=0; i<m*n; i++)
  {
    if (!CMP(truth[i], test[i]))
    {
      printf("FAIL\n");
      
      if (printonfail)
      {
        printf("Test matrix:\n");
        matprinter(m, n, test);
        printf("\nTruth matrix:\n");
        matprinter(m, n, truth);
      }
      return -1;
    }
  }
  
  printf("PASS\n");
  
  return 0;
}

#endif
