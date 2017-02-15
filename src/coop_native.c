/* Automatically generated. Do not edit by hand. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP R_check_badvals(SEXP x_);
extern SEXP R_co_mat(SEXP x, SEXP type_, SEXP inplace_, SEXP trans_, SEXP inv_);
extern SEXP R_co_mat_pairwise(SEXP x, SEXP type_, SEXP inv_);
extern SEXP R_co_matmat(SEXP x, SEXP y, SEXP type_, SEXP inplace_, SEXP trans_, SEXP inv_);
extern SEXP R_co_sparse(SEXP n_, SEXP a, SEXP i, SEXP j, SEXP index_, SEXP type_, SEXP inv_);
extern SEXP R_co_vecvec(SEXP x, SEXP y, SEXP type_);
extern SEXP R_csc_to_coo(SEXP row_ind, SEXP col_ptr);
extern SEXP R_fast_naomit(SEXP x);
extern SEXP R_fast_naomit_dbl(SEXP x_);
extern SEXP R_fast_naomit_int(SEXP x_);
extern SEXP R_naomit_vecvec(SEXP x_, SEXP y_);
extern SEXP R_scaler(SEXP centerx_, SEXP scalex_, SEXP x);
extern SEXP R_sparsity_dbl(SEXP x, SEXP tol);
extern SEXP R_sparsity_int(SEXP x);

static const R_CallMethodDef CallEntries[] = {
  {"R_check_badvals", (DL_FUNC) &R_check_badvals, 1},
  {"R_co_mat", (DL_FUNC) &R_co_mat, 5},
  {"R_co_mat_pairwise", (DL_FUNC) &R_co_mat_pairwise, 3},
  {"R_co_matmat", (DL_FUNC) &R_co_matmat, 6},
  {"R_co_sparse", (DL_FUNC) &R_co_sparse, 7},
  {"R_co_vecvec", (DL_FUNC) &R_co_vecvec, 3},
  {"R_csc_to_coo", (DL_FUNC) &R_csc_to_coo, 2},
  {"R_fast_naomit", (DL_FUNC) &R_fast_naomit, 1},
  {"R_fast_naomit_dbl", (DL_FUNC) &R_fast_naomit_dbl, 1},
  {"R_fast_naomit_int", (DL_FUNC) &R_fast_naomit_int, 1},
  {"R_naomit_vecvec", (DL_FUNC) &R_naomit_vecvec, 2},
  {"R_scaler", (DL_FUNC) &R_scaler, 3},
  {"R_sparsity_dbl", (DL_FUNC) &R_sparsity_dbl, 2},
  {"R_sparsity_int", (DL_FUNC) &R_sparsity_int, 1},
  {NULL, NULL, 0}
};
void R_init_coop(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
