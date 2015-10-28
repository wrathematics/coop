#' fastcosim
#'
#' @name fastcosim-package
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
#' cosine similarities of the columns are computed.
#' 
#' 
#' @useDynLib fastcosim, R_cosine_mat, R_cosine_vecvec
#' 
#' @docType package
#' @author Drew Schmidt and Wei-Chen Chen
#' @keywords package
NULL
