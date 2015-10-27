#' fastcosim
#'
#' @name fastcosim-package
#' 
#' @useDynLib fastcosim, cosine_fill_loop
#' 
#' @docType package
#' @author Drew Schmidt and Wei-Chen Chen
#' @keywords package
NULL



#' Cosine Similarity
#' 
#' Compute the cosine similarity matrix efficiently.
#' 
#' @details
#' The function computes the cosine similarity between all column
#' vectors of the numeric input matrix
#' 
#' @param x
#' A numeric matrix.
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
cosine <- function(x)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(cosine_fill_loop, x)
}

