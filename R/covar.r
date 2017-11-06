#' Covariance
#' 
#' An optimized, efficient implemntation for computing covariance.
#' 
#' @details
#' See \code{?coop-package} for implementation details.
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
#' @name covar
#' @rdname covar
NULL

#' @rdname covar
#' @export
covar <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE) UseMethod("covar")



#' @export
covar.matrix <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  co_matrix(x, y, CO_VAR, use, inplace, trans=FALSE, inverse=inverse)
}

#' @export
covar.data.frame <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  x_mat = as.matrix(x)
  if (!missing(y) && is.data.frame(y))
    y = as.matrix(y)
  
  covar.matrix(x_mat, y, use=use, inplace=inplace, inverse=inverse)
}



#' @export
covar.default <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  co_vecvec(x, y, CO_VAR, use)
}



### tcovar

#' @rdname covar
#' @export
tcovar <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE) UseMethod("tcovar")



#' @export
tcovar.matrix <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  co_matrix(x, y, CO_VAR, use, inplace, trans=TRUE, inverse=inverse)
}

#' @export
tcovar.data.frame <- function(x, y, use="everything", inplace=FALSE, inverse=FALSE)
{
  x_mat = as.matrix(x)
  if (!missing(y) && is.data.frame(y))
    y = as.matrix(y)
  
  tcovar.matrix(x_mat, y, use=use, inplace=inplace, inverse=inverse)
}
