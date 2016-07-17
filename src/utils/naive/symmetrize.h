#ifndef __COOP_NAIVE_SYMMETRIZE_H__
#define __COOP_NAIVE_SYMMETRIZE_H__

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

#endif
