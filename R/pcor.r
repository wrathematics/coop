#' Pearson Correlation
#' 
#' An optimized, efficient implemntation for computing the pearson
#' correlation.
#' 
#' @details
#' See \code{?coop} for implementation details.
#' 
#' @param x
#' A numeric dataframe/matrix or vector.
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
#' @param inverse
#' Logical; should the inverse covariance matrix be returned?
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
#' @name pcor
#' @rdname pcor
NULL

#' @rdname pcor
#' @export
pcor <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE) UseMethod("pcor")



#' @export
pcor.matrix <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  co_matrix(x, y, CO_ORR, use, inplace, trans=FALSE, inverse=inverse)
}

#' @export
pcor.data.frame <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  x_mat = as.matrix(x)
  if (!missing(y) && is.data.frame(y))
    y = as.matrix(y)
  
  pcor.matrix(x_mat, y, use=use, inplace=inplace, inverse=inverse)
}



#' @export
pcor.default <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  co_vecvec(x, y, CO_ORR, use)
}



### tpcor

#' @rdname pcor
#' @export
tpcor <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE) UseMethod("tpcor")



#' @export
tpcor.matrix <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  co_matrix(x, y, CO_ORR, use, inplace, trans=TRUE, inverse=inverse)
}

#' @export
tpcor.data.frame <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  x_mat = as.matrix(x)
  if (!missing(y) && is.data.frame(y))
    y = as.matrix(y)
  
  tpcor.matrix(x_mat, y, use=use, inplace=inplace, inverse=inverse)
}
