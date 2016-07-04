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

// Functions for computing covariance, (pearson) correlation, and cosine similarity

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "coop.h"
#include "utils/safeomp.h"


#if 0


#define WT_UNBIASED 1
#define WT_ML       2

#define BADWT -1

static inline wtchecks(const int m, const double *wt))
{
  double sum = 0;
  
  SAFE_FOR_SIMD
  for (int i=0; i<m; i++)
  {
    if (wt[i] < 0)
      return BADWT;
    
    sum += wt[i];
  }
  
  if (sum != 1.0)
    return BADWT;
}

static inline cp(const int m, const int n, const double alpha, const double * const restrict x, double * restrict c)
{
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &alpha, x, &m, &(double){0.0}, c, &n);
}

static void wtcp(const int m, const int n, double * const restrict x, const int wtlen, const double * const restrict wt)
{
  double alpha;
  if (method == WT_UNBIASED)
  {
    if (wtlen == 1)
      alpha = 1. / (1. - ((double)m)*wt[0]*wt[0]);
    else
    {
      alpha = 0.;
      SAFE_FOR_SIMD
      for (int i=0; i<m; i++)
        alpha = wt[i]*wt[i];
      
      alpha = 1. - alpha;
    }
  }
  else
    alpha = 1.;
  
  cp(m, n, alpha, x, c);
}



static void center_wt(const int m, const int n, const double * const restrict x, const double * const restrict wt, double * restrict colmeans)
{
  #pragma omp parallel for default(none) if(m*n>OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    const int mj = m*j;
    colmeans[j] = 0.;
    
    SAFE_FOR_SIMD
    for (int i=0; i<m; i++)
    {
      colmeans[j] += wt[j] * x[i + mj];
    }
  }
}



// TODO just operate on x in place, leave the copy to the user
int coop_covar_wt_mat(const int m, const int n, const double * const restrict x, int wtlen, const double * const restrict wt, double * restrict colmeans, double *restrict cov)
{
  double wtval;
  double wt_pt;
  
  if (wt == NULL)
  {
    wtval = 1./((double) m);
    wt_pt = &wtval;
    wtlen = 1;
  }
  else
    wt_pt = wt;
  
  
  double *x_cp = malloc(m*n*sizeof(*x));
  CHECKMALLOC(x_cp);
  memcpy(x_cp, x, m*n*sizeof(*x));
  
  center_wt(m, n, x, wt, colmeans);
  
  coop_scale(true, false, m, n, x_cp, colmeans, NULL);
  
  wtcp(m, n, x, wt)
  
  coop_symmetrize(n, cov);
  
  return 0;
}



n <- nrow(x)
x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
cov <- switch(match.arg(method), unbiased = crossprod(x)/(1 - 
    sum(wt^2)), ML = crossprod(x))

y <- list(cov = cov, center = center, n.obs = n)
if (with.wt) 
    y$wt <- wt
if (cor) {
    Is <- 1/sqrt(diag(cov))
    R <- cov
    R[] <- Is * cov * rep(Is, each = nrow(cov))
    y$cor <- R
}
y


#endif
