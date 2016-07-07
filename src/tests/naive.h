#ifndef __COOP_TESTS_NAIVE_H__
#define __COOP_TESTS_NAIVE_H__

static inline void matprinter(int m, int n, double *x)
{
  int i, j;
  
  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
      printf("%.0f ", x[i+m*j]);
    
    putchar('\n');
  }
  
  putchar('\n');
}

static inline void symmetrize_naive(const int n, double *restrict x)
{
  #pragma omp parallel for default(none) shared(x) schedule(dynamic, 1) if(n>OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    const int nj = n*j;
    
    for (int i=j+1; i<n; i++)
      x[j + n*i] = x[i + nj];
  }
}

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

#endif
