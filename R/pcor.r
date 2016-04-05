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
#' @param use
#' The NA handler, as in R's \code{cov()} and \code{cor()}
#' functions.  Options are "everything", "all.obs", and 
#' "complete.obs".
#' @param inplace
#' Logical; if \code{TRUE} then the method used is slower but
#' uses less memory than if \code{FALSE}.  See \code{?coop-package}
#' for details.
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
pcor <- function(x, y, use="everything", inplace=FALSE) UseMethod("pcor")



#' @export
pcor.matrix <- function(x, y, use="everything", inplace=FALSE)
{
  co_matrix(x, y, CO_ORR, use, inplace)
}



#' @export
pcor.default <- function(x, y, use="everything", inplace=FALSE)
{
  co_vecvec(x, y, CO_ORR, use)
}
