/*  Copyright (c) 2016, Schmidt
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <stdlib.h>
#include <string.h>

#include "coop.h"
#include "utils/inverse.h"
#include "utils/lapack.h"
#include "utils/mmult.h"
#include "utils/safeomp.h"
#include "utils/sumstats.h"


// matrix-matrix multiply: diagonal * general
// di is mxm, x is mxn
static inline void mm_diXge(const bool inv, const int m, const int n, const double *const restrict di, double *restrict x)
{
  if (inv)
  {
    #pragma omp parallel for
    for (int j=0; j<n; j++)
    {
      const register int mj = m*j;
      SAFE_SIMD
      for (int i=0; i<m; i++)
      {
        const double tmp = 1.0 / di[i];
        x[i + mj] *= tmp;
      }
    }
  }
  else
  {
    #pragma omp parallel for
    for (int j=0; j<n; j++)
    {
      const register int mj = m*j;
      SAFE_SIMD
      for (int i=0; i<m; i++)
      {
        const double tmp = di[i];
        x[i + mj] *= tmp;
      }
    }
  }
}




// ---------------------------------------------
//  QR Methods
// ---------------------------------------------


int cov_qr(const bool takeinv, const int m, const int n, double *restrict x, double *restrict cov)
{
  int info;
  int lwork;
  double tmp;
  const char uplo = 'u';
  const char trans = 't';
  const char diag = 'n';
  char side;
  double alpha;
  double *tau, *work, *x_cp;
  
  // allocations: tau, work, x_cp
  tau = malloc(MAX(m, n) * sizeof(*tau));
  if (tau == NULL)
    return COOP_BADMALLOC;
  
  lwork = -1;
  dgeqrf_(&m, &n, x, &m, tau, &tmp, &lwork, &info);
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  if (work == NULL)
  {
    info = COOP_BADMALLOC;
    goto allocfail_work;
  }
  
  x_cp = malloc(m*n*sizeof(*x));
  if (x_cp == NULL)
  {
    info = COOP_BADMALLOC;
    goto allocfail_xcp;
  }
  
  // center
  memcpy(x_cp, x, m*n*sizeof(*x));
  remove_colmeans(m, n, x_cp);
  
  // X = QR
  dgeqrf_(&m, &n, x, &m, tau, work, &lwork, &info);
  if (info != COOP_OK)
    goto cleanup;
  
  
  // compute (xtx)^-1 = x^-1 x^T^-1 or xtx
  memset(cov, 0.0, n*n * sizeof(*cov));
  dlacpy_(&uplo, &n, &n, x, &m, cov, &n);
  
  if (takeinv)
  {
    side = 'r';
    alpha = ((double) (m-1));
    inv_tri(m, n, x);
    // dtrmm_(&side, &uplo, &trans, &diag, &n, &n, &alpha, x, lda, b, ldb);
  }
  else
  {
    side = 'l';
    alpha = 1. / ((double) (m-1));
    // dtrmm_(&side, &uplo, &trans, &diag, &n, &n, &alpha, x, lda, b, ldb);
  }
  
  
  cleanup:
  free(x_cp);
  allocfail_xcp:
  free(work);
  allocfail_work:
  free(tau);
  
  return info;
}



// ---------------------------------------------
//  SVD Methods
// ---------------------------------------------


static int svd_nou(const int m, const int n, double *restrict x, double *restrict s, double *restrict vt)
{
  char jobz;
  int i;
  int info = 0;
  int lwork, *iwork;
  double tmp, *work;
  double *u;
  const int min_mn = m<n ? m : n;
  
  jobz = 's';
  u = malloc(min_mn*n * sizeof(*u));
  
  
  
  iwork = malloc(8*min_mn * sizeof(*iwork));
  
  lwork = -1;
  dgesdd_(&jobz, &m, &n, x, &m, s, u, &m, vt, &min_mn, &tmp, &lwork, iwork, &info);
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  dgesdd_(&jobz, &m, &n, x, &m, s, u, &m, vt, &min_mn, work, &lwork, iwork, &info);
  
  free(work);
  free(iwork);
  if (u != NULL)
    free(u);
  
  return info;
}



/*
int covinv_svd()
{
  // factor svd
  
  // sigma = (sigma^2)^{-1}
  inv_diagsq(n, s);
  
  // out = s*vt
  mm_diXge(n, n, s, vt, out);
  
  // out = v * out
  
}
*/



int cov_svd(const bool takeinv, const int m, const int n, double *restrict x, double *restrict cov)
{
  int info;
  int lwork, tmp;
  char side;
  double alpha;
  double *s, *vt, *work, *x_cp;
  char transx, transy;
  
  const int min_mn = (m<n?m:n);
  
  // allocations: 
  s = malloc(min_mn * sizeof(*s));
  vt = malloc(min_mn*n * sizeof(*vt));
  
  x_cp = malloc(m*n*sizeof(*x));
  if (x_cp == NULL)
  {
    info = COOP_BADMALLOC;
    goto allocfail_xcp;
  }
  
  // center
  memcpy(x_cp, x, m*n*sizeof(*x));
  remove_colmeans(m, n, x_cp);
  
  // take svd
  info = svd_nou(m, n, x_cp, s, vt);
  
  if (info != COOP_OK)
    goto cleanup;
  
  mm_diXge(takeinv, m, n, s, vt);
  
  matmult(false, true, 1.0, min_mn, n, vt, min_mn, n, vt, cov);
  
  cleanup:
  free(x_cp);
  allocfail_xcp:
  
  return info;
}
