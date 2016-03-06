#' Covariance
#' 
#' Compute covariance efficiently.
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
#' The covariance.
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
covar <- function(x, y) UseMethod("covar")



#' @export
covar.matrix <- function(x, y)
{
  co_matrix(x, y, CO_VAR)
}



#' @export
covar.default <- function(x, y)
{
  co_vecvec(x, y, CO_VAR)
}
