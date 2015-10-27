#' fastcosim
#'
#' @name fastcosim-package
#' 
#' @useDynLib fastcosim, R_cosine_mat, R_cosine_vecvec
#' 
#' @docType package
#' @author Drew Schmidt and Wei-Chen Chen
#' @keywords package
NULL



#' Cosine Similarity
#' 
#' Compute the cosine similarity matrix efficiently.  The function
#' syntax and behavior is largely modeled after that of the
#' \code{cosine()} function from the \code{lsa} package, although
#' with a very different implementation.
#' 
#' @details
#' The function computes the cosine similarity between all column
#' vectors of the numeric input matrix
#' 
#' @param x
#' A numeric matrix or vector.
#' @param y
#' A vector (when \code{x} is a vector) or missing (blank) when 
#' \code{x} is a matrix.
#' 
#' @return
#' The \eqn{n\times n} matrix of all pair-wise vector cosine
#' similarities of the columns.  The implementation is dominated
#' by a symmetric rank-k update via the BLAS function dsyrk(); see
#' the package \code{README.md} file for a discussion of the 
#' algorithm implementation and complexity.
#' 
#' @examples
#' library(fastcosim)
#' x <- matrix(rnorm(10*3), 10, 3)
#' cosine(x)
#' 
#' @export
cosine <- function(x, y)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (missing(y))
  {
    .Call(R_cosine_mat, x)
  }
  else
  {
    if (!is.double(y))
      storage.mode(y) <- "double"
    
    .Call(R_cosine_vecvec, x, y)
  }
}

