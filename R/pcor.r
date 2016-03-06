#' Pearson Correlation
#' 
#' Compute the pearson correlation efficiently.
#' 
#' @details
#' This implementation uses the internals for \code{cosine()}.
#' 
#' @param x
#' A numeric matrix or vector.
#' @param y
#' A vector (when \code{x} is a vector) or missing (blank) when 
#' \code{x} is a matrix.
#' 
#' @return
#' The \eqn{n\times n} matrix of all pair-wise vector pearson
#' correlations of the columns.
#' 
#' @examples
#' library(fastco)
#' x <- matrix(rnorm(10*3), 10, 3)
#' pcor(x)
#' 
#' pcor(x[, 1], x[, 2])
#' 
#' @seealso \code{\link{cosine}}
#' @export
pcor <- function(x, y) UseMethod("pcor")



#' @export
pcor.matrix <- function(x, y)
{
  co_matrix(x, y, CO_ORR)
}



#' @export
pcor.default <- function(x, y)
{
  co_vecvec(x, y, CO_ORR)
}
