#' Covariance
#' 
#' An optimized, efficient implemntation for computing covariance.
#' 
#' @details
#' See \code{?coop-package} for implementation details.
#' 
#' @param x
#' A numeric matrix or vector.
#' @param y
#' A vector (when \code{x} is a vector) or missing (blank) when 
#' \code{x} is a matrix.
#' 
#' @return
#' The covariance matrix.
#' 
#' @examples
#' x <- matrix(rnorm(10*3), 10, 3)
#' 
#' coop::pcor(x)
#' coop::pcor(x[, 1], x[, 2])
#' 
#' @author Drew Schmidt
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
