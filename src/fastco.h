/*  Copyright (c) 2015-2016, Schmidt
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

#ifndef FASTCOSIM_H
#define FASTCOSIM_H


#define TMP_VEC_SIZE 1024

#define EPSILON 1e-10


#define CHECKMALLOC(x) if(x==NULL) return -1

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

// BLAS
void dgemm_(const char *transa, const char *transb, const int *m, const int *n, 
            const int *k, const double *restrict alpha, const double *restrict a, 
            const int *lda, const double *restrict b, const int *ldb, 
            const double *beta, double *restrict c, const int *ldc);

void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k, 
            const double *restrict alpha, const double *restrict a, const int *lda, 
            const double *restrict beta, double *restrict c, const int *ldc);

// dense
void cosine_mat(const int m, const int n, const double *restrict x, double *restrict cos);
double cosine_vecvec(const int n, const double *restrict x, const double *restrict y);

void pcor_mat(const int m, const int n, const double *restrict x, double *restrict cor);
double pcor_vecvec(const int n, const double *restrict x, const double *restrict y);

void covar_mat(const int m, const int n, const double *restrict x, double *restrict cov);
double covar_vecvec(const int n, const double *x, const double *y);

// sparse
int cosine_sparse_coo(const int index, const int n, const int len, const double *restrict a, const int *restrict rows, const int *restrict cols, double *restrict cos);

// utils
void diag2one(const unsigned int n, double *restrict x);
void symmetrize(const int n, double *restrict x);


#endif
