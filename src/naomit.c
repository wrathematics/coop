#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

#include "omp.h"


static SEXP copymat(const int m, const int n, SEXP x_)
{
  SEXP ret;
  int i;
  
  if (m*n < OMP_MIN_SIZE)
  {
    ret = PROTECT(duplicate(x_));
  }
  else
  {
    const double *x = REAL(x_);
    
    PROTECT(ret = allocMatrix(REALSXP, m, n));
    double *retptr = REAL(ret);
    
    #pragma omp parallel for simd default(shared) if(m*n>OMP_MIN_SIZE)
    for (i=0; i<m*n; i++)
      retptr[i] = x[i];
  }
  
  UNPROTECT(1);
  return ret;
}



SEXP R_fast_naomit_dbl(SEXP x_)
{
  int i, j, mj;
  const int m = nrows(x_);
  const int n = ncols(x_);
  SEXP ret;
  int *rows = (int*) calloc(m, sizeof(*rows));
  int m_fin = m;
  int row;
  int itmp;
  
  const double *x = REAL(x_);
  
  #pragma omp parallel for default(shared) private(i, j) if(m*n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    mj = m*j;
    
    SAFE_SIMD
    for (i=0; i<m; i++)
    {
      if (ISNA(x[i + mj]) || ISNAN(x[i + mj]))
        rows[i] = 1;
    }
  }
  
  SAFE_FOR_SIMD
  for (i=0; i<m; i++)
    m_fin -= rows[i];
  
  if (m_fin == m)
  {
    ret = copymat(m, n, x_);
    free(rows);
    return ret;
  }
  
  PROTECT(ret = allocMatrix(REALSXP, m_fin, n));
  double *retptr = REAL(ret);
  
  #pragma omp parallel for default(shared) private(i, j, row, itmp, mj) if(m*n>OMP_MIN_SIZE)
  for (j=0; j<n; j++)
  {
    mj = m*j;
    row = 0;
    
    SAFE_SIMD
    for (i=0; i<m; i++)
    {
      itmp = rows[i];
      if (!itmp)
      {
        retptr[row + m_fin*j] = x[i + mj];
        row++;
      }
    }
  }
  
  free(rows);
  UNPROTECT(1);
  return ret;
}
