#' fastcosim
#'
#' @name fastcosim-package
#' 
#' @useDynLib fastcosim, cosine_fill_loop
#' 
#' @docType package
#' @author Drew Schmidt and Wei-Chen Chen
#' @keywords package
NULL



#' Cosine Similarity
#' 
#' Compute cosine similarity efficiently.
#' 
#' @param x
#' A numeric matrix.
#' 
#' @export
cosine <- function(x)
{
  cp <- crossprod(x)
  dg <- diag(cp)
  
  .Call(cosine_fill_loop, cp, dg)
}

