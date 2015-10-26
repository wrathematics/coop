#include <R.h>
#include <Rinternals.h>
#include <math.h>


static void fill(const unsigned int n, double *restrict cosim, double *restrict crossprod, double *restrict diag)
{
  int i, j;
  double diagj;
  
  // Fill lower triangle and diagonal
  for (j=0; j<n; j++)
  {
    diagj = diag[j];
    
    for (i=0; i<j; i++)
      cosim[i + n*j] = crossprod[i + n*j] / sqrt(diag[i] * diagj);
    
    cosim[i + n*i] = 1.0;
  }
  
  // Copy lower triangle to upper
  for (j=0; j<n; j++)
  {
    for (i=j+1; i<n; i++)
      cosim[i + n*j] = cosim[j + n*i];
  }
}



SEXP cosine_fill_loop(SEXP crossprod, SEXP diag)
{
  const unsigned int n = nrows(crossprod);
  SEXP cosim;
  PROTECT(cosim = allocMatrix(REALSXP, n, n));
  double *cosim_pt = REAL(cosim);
  double *crossprod_pt = REAL(crossprod);
  double *diag_pt = REAL(diag);
  
  fill(n, cosim_pt, crossprod_pt, diag_pt);
  
  UNPROTECT(1);
  return cosim;
}

