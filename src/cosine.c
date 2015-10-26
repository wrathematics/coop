#include <R.h>
#include <Rinternals.h>
#include <math.h>

SEXP cosine_fill_loop(SEXP crossprod, SEXP diag)
{
  const unsigned int n = nrows(crossprod);
  SEXP cosim;
  PROTECT(cosim = allocMatrix(REALSXP, n, n));
  
  // Fill lower triangle and diagonal
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<=j; i++)
    {
      if (i == j)
        REAL(cosim)[j + n*j] = 1.0;
      else
        REAL(cosim)[i + n*j] = REAL(crossprod)[i + n*j] / sqrt(REAL(diag)[i] * REAL(diag)[j]);
    }
  }
  
  // Copy lower triangle to upper
  for (int j=0; j<n; j++)
  {
    for (int i=j+1; i<n; i++)
      REAL(cosim)[i + n*j] = REAL(cosim)[j + n*i];
  }
  
  UNPROTECT(1);
  return cosim;
}

