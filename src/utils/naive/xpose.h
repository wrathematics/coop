#ifndef __COOP_NAIVE_XPOSE_H__
#define __COOP_NAIVE_XPOSE_H__

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
