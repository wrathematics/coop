#' fastco
#'
#' @name fastco-package
#' 
#' @description
#' A micro-package for computing cosine similarity.  This is
#' commonly needed in, for example, natural language processing,
#' where the cosine similarity coefficients of all columns of a
#' term-document or document-term matrix is needed.
#' 
#' @section Implementation Details:
#' The only exported function of this package, \code{cosine()},
#' takes two forms.  If supplied with two vectors, then the cosine
#' similarity of the two vectors will be computed.  On the other
#' hand, if supplied a single matrix, then all pair-wise vector
#' cosine similarities of the columns are computed.  The storage
#' options for a "matrix" are ordinary R matrices and "simple
#' triplet matrix" from the slam package.
#' 
#' 
#' @useDynLib fastco, R_cosine_mat, R_cosine_vecvec,
#'   R_pcor_mat, R_pcor_vecvec,
#'   R_cosine_sparse_coo, R_sparsity_int, R_sparsity_dbl
#' 
#' @docType package
#' @author Drew Schmidt and Wei-Chen Chen
#' @keywords package
NULL
