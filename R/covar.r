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
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!missing(y))
    stop("argument 'y' can not be used with a matrix 'x'")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(R_co_mat, x, 3L)
}



#' @export
covar.default <- function(x, y)
{
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  
  if (missing(y))
    return(1.0)
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  
  if (length(x) != length(y))
    stop("vectors 'x' and 'y' must have the same length")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  .Call(R_co_vecvec, x, y, 3L)
}
