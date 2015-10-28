#' fastcosim
#'
#' @name fastcosim-package
#' 
#' @description
#' A micro-package for computing cosine similarity.  This is
#' commonly needed in, for example, natural language processing,
#' where the cosine similarity coefficients of all columns of a
#' term-document or document-term matrix is needed.
#' 
#' @section Implementation Details:
#' The only exported function of this package, \code{cosine()},
#' takes two forms.  If supplied with two vectors, then the cosine
#' similarity of the two vectors will be computed.  On the other
#' hand, if supplied a single matrix, then all pair-wise vector
#' cosine similarities of the columns are computed.
#' 
#' 
#' @useDynLib fastcosim, R_cosine_mat, R_cosine_vecvec
#' 
#' @docType package
#' @author Drew Schmidt and Wei-Chen Chen
#' @keywords package
NULL



#' Cosine Similarity
#' 
#' Compute the cosine similarity matrix efficiently.  The function
#' syntax and behavior is largely modeled after that of the
#' \code{cosine()} function from the \code{lsa} package, although
#' with a very different implementation.
#' 
#' @details
#' The function computes the cosine similarity between all column
#' vectors of the numeric input matrix, or of two vectors, depending
#' on the input argument(s).
#' 
#' The implementation for matrix inputs is dominated
#' by a symmetric rank-k update via the BLAS subroutine \code{dsyrk};
#' see the package \code{README.md} file for a discussion of the 
#' algorithm implementation and complexity.
#' 
#' The implementation for two vector inputs is the product
#' \code{t(x) \%*\% y} via the BLAS subroutine \code{dgemm} divided
#' by the square root of the product of the crossproducts
#' \code{t(x) \%*\% x} and \code{t(y) \%*\% y} each computed via the
#' BLAS function \code{dsyrk()}.
#' 
#' @param x
#' A numeric matrix or vector.
#' @param y
#' A vector (when \code{x} is a vector) or missing (blank) when 
#' \code{x} is a matrix.
#' 
#' @return
#' The \eqn{n\times n} matrix of all pair-wise vector cosine
#' similarities of the columns.
#' 
#' @examples
#' library(fastcosim)
#' x <- matrix(rnorm(10*3), 10, 3)
#' cosine(x)
#' 
#' cosine(x[, 1], x[, 2])
#' 
#' @export
cosine <- function(x, y)
{
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  
  if (missing(y))
  {
    if (!is.matrix(x))
      stop("argument 'x' must be a matrix")
    
    if (!is.double(x))
      storage.mode(x) <- "double"
    
    .Call(R_cosine_mat, x)
  }
  else
  {
    if (!is.numeric(y))
      stop("argument 'y' must be numeric")
    
    if (length(x) != length(y))
      stop("vectors 'x' and 'y' must have the same length")
    
    if (!is.double(x))
      storage.mode(x) <- "double"
    if (!is.double(y))
      storage.mode(y) <- "double"
    
    .Call(R_cosine_vecvec, x, y)
  }
}

