#' Weighted Co-Operation
#' 
#' An optimized, efficient implemntation for computing weighted covariance,
#' correlation, and cosine similarity.  Similar to R's \code{cov.wt()}.
#' 
#' @details
#' See \code{?coop-package} for implementation details.
#' 
#' @param x
#' A matrix or data.frame.
#' @param wt
#' 
#' @param method
#' Either "unbiased" or "ml". Unlike R, case is ignored.
#' 
#' @examples
#' x <- matrix(rnorm(10*3), 10, 3)
#' cov.wt(x)
#' 
#' @author Drew Schmidt
#' @seealso \code{\link{cosine}}, \code{\link{pcor}}, and \code{\link{covar}}
#' @name weighted
#' @rdname weighted
NULL

#' @rdname weighted
#' @export
cosine_wt <- function(x, wt=NULL, method="unbiased") UseMethod("cosine_wt")

#' @rdname weighted
#' @export
pcor_wt <- function(x, wt=NULL, method="unbiased") UseMethod("pcor_wt")

#' @rdname weighted
#' @export
covar_wt <- function(x, wt=NULL, method="unbiased") UseMethod("covar_wt")




co_wt.data.frame <- function(x, wt=NULL, method="unbiased")
{
  covar_wt.matrix(as.matrix(x), wt=wt, method=method)
}

cosine_wt.data.frame <- co_wt.data.frame
pcor_wt.data.frame <- co_wt.data.frame
covar_wt.data.frame <- co_wt.data.frame



cosine_wt.matrix <- function(x, wt=NULL, method="unbiased")
{
  co_wt(x=x, wt=wt, method=method, type=CO_SIM)
}



pcor_wt.matrix <- function(x, wt=NULL, method="unbiased")
{
  co_wt(x=x, wt=wt, method=method, type=CO_ORR)
}



covar_wt.matrix <- function(x, wt=NULL, method="unbiased")
{
  covar_wt(x=x, wt=wt, method=method, type=CO_VAR)
}



co_wt <- function(x, wt=NULL, method="unbiased", type)
{
  method <- match.arg(tolower(method), c("unbiased", "ml"))
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  check_badvals(x)
  
  .Call(R_cov_wt, x, wt, type=type)
}
