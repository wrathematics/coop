#' scaler
#' 
#' A function to center (subtract mean) and/or scale (divide by
#' standard deviation) data column-wise in a computationally
#' efficient way.
#' 
#' @details
#' Unlike its R counterpart \code{scale()}, the arguments
#' \code{center} and \code{scale} can only be logical values
#' (and not vectors).
#' 
#' @param x
#' The input matrix.
#' @param center,scale
#' Logical; determine if the data should be centered and/or scaled.
#' 
#' @return
#' The centered/scaled data, with attributes as in R's \code{scale()}.
#' 
#' @export
scaler <- function(x, center=TRUE, scale=TRUE)
{
  check.is.flag(center)
  check.is.flag(scale)
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!center && !scale)
    return(x)
  else
    .Call(R_scaler, center, scale, x)
}
