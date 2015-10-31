#' Sparsity
#' 
#' Show the sparsity (as a count or proportion) of a matrix.  For
#' example, .99 sparsity means 99\% of the values are zero.
#' Similarly, a sparsity of 0 means the matrix is fully dense.
#' 
#' @details
#' The implementation is very efficient for dense matrices.  For
#' sparse triplet matrices, the count is trivial.
#' 
#' @param x
#' The matrix, stored as an ordinary R matrix or as a "simple
#' triplet matrix" (from the slam package).
#' @param proportion
#' Logical; determines if a proportion or a count should be returned.
#' 
#' @return
#' The sparsity of the input matrix, of the form a proportion or
#' a count.
#' 
#' @examples
#' ## Completely sparse matrix
#' x <- matrix(0, 10, 10)
#' sparsity(x)
#' 
#' ## 15\% density / 85\% sparsity
#' x[sample(length(x), size=15)] <- 1
#' sparsity(x)
#' 
#' @export
sparsity <- function(x, proportion=TRUE) UseMethod("sparsity")



#' @export
sparsity.matrix <- function(x, proportion=TRUE)
{
  if (is.integer(x))
    count <- .Call(R_sparsity_int, x)
  else if (is.double(x))
    count <- .Call(R_sparsity_dbl, x, tol=1e-10)
  else
    stop("matrix 'x' must be numeric.")
  
  if (proportion)
    count / nrow(x) / ncol(x)
  else
    count
}



#' @export
sparsity.simple_triplet_matrix <- function(x, proportion=TRUE)
{
  count <- x$nrow*x$ncol - length(x$v)
  if (proportion)
    count / nrow(x) / ncol(x)
  else
    count
}

