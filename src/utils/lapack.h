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

#ifndef __COOP_LIB_LAPACK_H__
#define __COOP_LIB_LAPACK_H__

// BLAS

void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
            const int *k, const double *restrict alpha, const double *restrict a,
            const int *lda, const double *restrict b, const int *ldb,
            const double *beta, double *restrict c, const int *ldc);

void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k,
            const double *restrict alpha, const double *restrict a, const int *lda,
            const double *restrict beta, double *restrict c, const int *ldc);


// LAPACK
void dgeqrf_(const int *m, const int *n, double *restrict a, const int *lda, 
             double *restrict tau, double *restrict work, int *lwork, int *info);

void dgesdd_(const char *jobz, const int *m, const int *n, double *a, 
             const int *lda, double *restrict s, double *restrict u, const int *restrict ldu, double *restrict vt, 
             const int *ldvt, double *restrict work, const int *lwork, int *iwork, int *info);

void dgetrf_(const int *m, const int *n, double *restrict a, const int *lda, 
             int *restrict ipiv, int *info);

void dgetri_(const int *n, double *restrict a, const int *lda, 
             int *restrict ipiv, double *work, int *lwork, int *info);

void dtrtri_(const char *uplo, const char *diag, const int *n, double *restrict a,
             const int *lda, const int *info);

void dlacpy_(const char *uplo, const int *m, const int *n, const double *restrict x, 
             const int *lda, double *restrict y, const int *ldb);

void dpotrf_(const char *uplo, const int *n, double *restrict a, const int *lda, int *info);

void dpotri_(const char *uplo, const int *n, double *restrict a, const int *lda, int *info);


#endif
