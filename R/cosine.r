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
#' on the input argument(s). Two different storage schemes are 
#' accepted for the matrix version.  For dense matrices, an ordinary
#' R matrix input is accepted.  For sparse matrices, a matrix in
#' COO format, namely \code{simple_triplet_matrix} from the slam
#' package, is accepted.
#' 
#' The implementation for dense matrix inputs is dominated
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
#' library(fastco)
#' x <- matrix(rnorm(10*3), 10, 3)
#' cosine(x)
#' 
#' cosine(x[, 1], x[, 2])
#' 
#' @seealso \code{\link{sparsity}}
#' @export
cosine <- function(x, y) UseMethod("cosine")



#' @export
cosine.matrix <- function(x, y)
{
  co_matrix(x, y, CO_SIM)
}



#' @export
cosine.default <- function(x, y)
{
  co_vecvec(x, y, CO_SIM)
}



#' @export
cosine.simple_triplet_matrix <- function(x, y)
{
  co_sparse(x, y, type=CO_SIM)
}
