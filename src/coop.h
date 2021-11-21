/*  Copyright (c) 2015-2016 Drew Schmidt
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

#ifndef __COOP_LIB_H__
#define __COOP_LIB_H__

#include <stdbool.h>

#include "utils/cdefs.h"


// dense
int coop_cosine_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cos);
int coop_tcosine_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cos);
int coop_cosine_matmat(const bool trans, const bool inv, const int m, const int n, const double * const restrict x, const double * const restrict y, double *restrict cos);
int coop_cosine_vecvec(const int n, const double * const restrict x, const double * const restrict y, double *restrict cos);

int coop_pcor_mat(const bool inv, const int m, const int n,
  const double * const restrict x, double *restrict cor);
int coop_tpcor_mat(const bool inv, const int m, const int n,
  const double * const restrict x, double *restrict cor);
int coop_pcor_matmat(const bool trans, const bool inv, const int m, const int n, const double * const restrict x, const double * const restrict y, double *restrict cor);
int coop_pcor_vecvec(const int n, const double * const restrict x, const double * const restrict y, double *restrict cor);

int coop_covar_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cov);
int coop_tcovar_mat(const bool inv, const int m, const int n,
  const double *const restrict x, double *restrict cov);
int coop_covar_matmat(const bool trans, const bool inv, const int m, const int n, const double * const restrict x, const double * const restrict y, double *restrict cov);
int coop_covar_vecvec(const int n, const double * const x, const double * const y, double *restrict cov);

// dense - inplace
int coop_covar_mat_inplace(const bool inv, const int m, const int n, const double *restrict x, double *restrict cov);
int coop_pcor_mat_inplace(const bool inv, const int m, const int n, const double * const restrict x, double *restrict cor);

// dense - pairwise complete observations (inplace)
int coop_cosine_mat_inplace_pairwise(const bool inv, const int m, const int n, const double * const restrict x, double *restrict cos);
int coop_pcor_mat_inplace_pairwise(const bool inv, const int m, const int n, const double * const restrict x, double *restrict cor);
int coop_covar_mat_inplace_pairwise(const bool inv, const int m, const int n, const double * const restrict x, double *restrict cov);

// scale
int coop_scale(const bool centerx, const bool scalex, const int m, const int n, double *restrict x, double *restrict colmeans, double *restrict colvars);

// sparse
int coop_cosine_sparse_coo(const bool inv, const int index, const int n, const int len, const double * const restrict a, const int *restrict rows, const int *restrict cols, double *restrict cos);


#endif
