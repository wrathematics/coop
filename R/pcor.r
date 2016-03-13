#' Pearson Correlation
#' 
#' An optimized, efficient implemntation for computing the pearson
#' correlation.
#' 
#' @details
#' See \code{?coop} for implementation details.
#' 
#' @param x
#' A numeric matrix or vector.
#' @param y
#' A vector (when \code{x} is a vector) or missing (blank) when 
#' \code{x} is a matrix.
#' 
#' @return
#' The pearson correlation matrix.
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
