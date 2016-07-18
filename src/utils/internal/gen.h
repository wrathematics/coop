#ifndef __COOP_UTILS_INTERNAL_GEN_H__
#define __COOP_UTILS_INTERNAL_GEN_H__

#include <stdlib.h>
#include <string.h>

#include "../cdefs.h"

#define EXITCHECK(check) if (check == BADMALLOC){fflush(stdout);fprintf(stderr, "--- ERROR: bad malloc\n");exit(BADMALLOC);}


static inline int zeromat(const int n, double **x)
{
  *x = malloc(n * sizeof(**x));
  CHECKMALLOC(x);
  
  memset(*x, 0.0, n*sizeof(**x));
  
  return 0;
}

static inline double* zeromat2(const int n)
{
  double *x;
  int ret = zeromat(n, &x);
  EXITCHECK(ret);
  
  return x;
}



static inline int cpalloc(const int n, const double *const restrict src, double **dest)
{
  *dest = malloc(n * sizeof(**dest));
  CHECKMALLOC(dest);
  
  memcpy(*dest, src, n*sizeof(**dest));
  
  return 0;
}

static inline double* cpalloc2(const int n, const double *const restrict src)
{
  double *dest;
  int ret = cpalloc(n, src, &dest);
  EXITCHECK(ret);
  
  return dest;
}



static inline int gen_boring(const int n, double **x)
{
  *x = malloc(n * sizeof(**x));
  CHECKMALLOC(x);
  
  for (int i=0; i<n; i++)
    (*x)[i] = (double) i+1;
  
  return 0;
}

static inline double* gen_boring2(const int n)
{
  double *x;
  int ret = gen_boring(n, &x);
  EXITCHECK(ret);
  
  return x;
}



static inline int gen_runif(const int n, double **x)
{
  *x = malloc(n * sizeof(**x));
  CHECKMALLOC(x);
  
  for (int i=0; i<n; i++)
    (*x)[i] = (double) rand() / RAND_MAX;
  
  return 0;
}

static inline double* gen_runif2(const int n)
{
  double *x;
  int ret = gen_runif(n, &x);
  EXITCHECK(ret);
  
  return x;
}


#endif
