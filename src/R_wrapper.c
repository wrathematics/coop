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

#include <R.h>
#include <Rinternals.h>
#include <math.h>

#include "coop.h"
#include "utils/inverse.h"
#include "utils/sparsity.h"


#define CO_SIM 1
#define CO_ORR 2
#define CO_VAR 3

#define BADTYPE() error("Invalid 'type' argument; please report this to the package author")

#define INT(x) INTEGER(x)[0]

static inline void checkret(const int check)
{
  if (check)
  {
    if (check == COOP_BADMALLOC)
      error("unable to allocate necessary memory");
    else
      error("Unable to compute inverse");
  }
}

// ---------------------------------------------
//  Dense
// ---------------------------------------------


SEXP R_co_mat(SEXP x, SEXP type_, SEXP inplace_, SEXP trans_, SEXP inv_)
{
  SEXP ret;
  int check;
  const int type = INT(type_);
  const unsigned int m = nrows(x);
  const unsigned int n = ncols(x);
  const int inplace = INT(inplace_);
  const int trans = INT(trans_);
  const int inv = INT(inv_);
  
  
  if (trans)
    PROTECT(ret = allocMatrix(REALSXP, m, m));
  else
    PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  
  if (type == CO_SIM)
    check = coop_cosine_mat(trans, inv, m, n, REAL(x), REAL(ret));
  else if (type == CO_ORR)
  {
    if (inplace)
      check = coop_pcor_mat_inplace(inv, m, n, REAL(x), REAL(ret));
    else
      check = coop_pcor_mat(trans, inv, m, n, REAL(x), REAL(ret));
  }
  else if (type == CO_VAR)
  {
    if (inplace)
      check = coop_covar_mat_inplace(inv, m, n, REAL(x), REAL(ret));
    else
      check = coop_covar_mat(trans, inv, m, n, REAL(x), REAL(ret));
  }
  else
    BADTYPE();
  
  UNPROTECT(1);
  
  checkret(check);
  
  return ret;
}



SEXP R_co_matmat(SEXP x, SEXP y, SEXP type_, SEXP inplace_, SEXP trans_, SEXP inv_)
{
  SEXP ret;
  int check;
  const int type = INT(type_);
  const unsigned int m = nrows(x);
  const unsigned int n = ncols(x);
  // const int inplace = INT(inplace_);
  const int trans = INT(trans_);
  const int inv = INT(inv_);
  
  
  if (trans)
    PROTECT(ret = allocMatrix(REALSXP, m, m));
  else
    PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  
  if (type == CO_SIM)
    check = coop_cosine_matmat(trans, inv, m, n, REAL(x), REAL(y), REAL(ret));
  else if (type == CO_ORR)
  {
    // if (inplace)
    // 
    // else
      check = coop_pcor_matmat(trans, inv, m, n, REAL(x), REAL(y), REAL(ret));
  }
  else if (type == CO_VAR)
  {
    // if (inplace)
    // 
    // else
      check = coop_covar_matmat(trans, inv, m, n, REAL(x), REAL(y), REAL(ret));
  }
  else
    BADTYPE();
  
  UNPROTECT(1);
  
  checkret(check);
  
  return ret;
}



SEXP R_co_vecvec(SEXP x, SEXP y, SEXP type_)
{
  SEXP ret;
  int check;
  const int type = INT(type_);
  const unsigned int n = LENGTH(x);
  double co;
  PROTECT(ret = allocVector(REALSXP, 1));
  
  if (type == CO_SIM)
    check = coop_cosine_vecvec(n, REAL(x), REAL(y), &co);
  else if (type == CO_ORR)
    check = coop_pcor_vecvec(n, REAL(x), REAL(y), &co);
  else if (type == CO_VAR)
    check = coop_covar_vecvec(n, REAL(x), REAL(y), &co);
  else
    BADTYPE();
  
  REAL(ret)[0] = co;
  UNPROTECT(1);
  
  checkret(check);
  
  return ret;
}



SEXP R_co_mat_pairwise(SEXP x, SEXP type_, SEXP inv_)
{
  SEXP ret;
  int check;
  const int type = INT(type_);
  const unsigned int m = nrows(x);
  const unsigned int n = ncols(x);
  const int inv = INT(inv_);
  PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  
  if (type == CO_SIM)
    check = coop_cosine_mat_inplace_pairwise(inv, m, n, REAL(x), REAL(ret));
  else if (type == CO_ORR)
    check = coop_pcor_mat_inplace_pairwise(inv, m, n, REAL(x), REAL(ret));
  else if (type == CO_VAR)
    check = coop_covar_mat_inplace_pairwise(inv, m, n, REAL(x), REAL(ret));
  else
    BADTYPE();
  
  UNPROTECT(1);
  
  checkret(check);
  
  return ret;
}



SEXP R_scaler(SEXP centerx_, SEXP scalex_, SEXP x)
{
  SEXP ret, cm, cv;
  const int m = nrows(x);
  const int n = ncols(x);
  const int centerx = INT(centerx_);
  const int scalex = INT(scalex_);
  int ptct;
  double *colmeans, *colvars;
  
  PROTECT(ret = allocMatrix(REALSXP, m, n));
  ptct = 1;
  memcpy(REAL(ret), REAL(x), m*n*sizeof(double));
  
  if (centerx)
  {
    PROTECT(cm = allocVector(REALSXP, n));
    colmeans = REAL(cm);
    ptct++;
  }
  else
  {
    cm = R_NilValue; // avoids a bogus compiler warning
    colmeans = NULL;
  }
  
  if (scalex)
  {
    PROTECT(cv = allocVector(REALSXP, n));
    colvars = REAL(cv);
    ptct++;
  }
  else
  {
    cv = R_NilValue; // same
    colvars = NULL;
  }
  
  coop_scale(centerx, scalex, m, n, REAL(ret), colmeans, colvars);
  
  if (centerx)
    setAttrib(ret, install("scaled:center"), cm);
  if (scalex)
    setAttrib(ret, install("scaled:scale"), cv);
  
  UNPROTECT(ptct);
  return ret;
}



// ---------------------------------------------
//  Sparse
// ---------------------------------------------

// #define INEDEX_FROM_1 1

SEXP R_co_sparse(SEXP n_, SEXP a, SEXP i, SEXP j, SEXP index_, SEXP type_, SEXP inv_)
{
  int check;
  const int n = INT(n_);
  const int index = INT(index_);
  const int type = INT(type_);
  const int inv = INT(inv_);
  SEXP ret;
  PROTECT(ret = allocMatrix(REALSXP, n, n));
  
  if (type == CO_SIM)
    check = coop_cosine_sparse_coo(inv, index, n, LENGTH(a), REAL(a), INTEGER(i), INTEGER(j), REAL(ret));
  else
    BADTYPE();
  
  UNPROTECT(1);
  
  checkret(check);
  
  return ret;
}





// ---------------------------------------------
//  Sparse utils
// ---------------------------------------------

SEXP R_sparsity_int(SEXP x)
{
  const int m = nrows(x), n = ncols(x);
  SEXP ret;
  
  PROTECT(ret = allocVector(INTSXP, 1));
  INT(ret) = coop_sparsity_int(m, n, INTEGER(x));
  UNPROTECT(1);
  
  return ret;
}



SEXP R_sparsity_dbl(SEXP x, SEXP tol)
{
  const int m = nrows(x), n = ncols(x);
  SEXP ret;
  
  PROTECT(ret = allocVector(INTSXP, 1));
  INT(ret) = coop_sparsity_dbl(m , n, REAL(x), REAL(tol)[0]);
  UNPROTECT(1);
  
  return ret;
}



#define INTi(x,i) INTEGER(x)[i]

SEXP R_csc_to_coo(SEXP row_ind, SEXP col_ptr)
{
  int j = 0;
  int c = 0, ind = 0;
  int diff;
  const int len = LENGTH(row_ind);
  
  SEXP col_ind;
  PROTECT(col_ind = allocVector(INTSXP, len));
  
  
  for (c=0; c<LENGTH(col_ptr)-1; c++) // hehehe
  {
    diff = INTi(col_ptr, c+1) - INTi(col_ptr, c);
    
    if (diff < 0)
      error("malformed dgCMatrix; impossible col_ptr array");
    
    while (diff > 0)
    {
      INTi(col_ind, ind) = j;
      
      ind++;
      diff--;
    }
    
    j++;
  }
  
  UNPROTECT(1);
  return col_ind;
}
