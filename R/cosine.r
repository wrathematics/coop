#' Cosine Similarity
#' 
#' Compute the cosine similarity matrix efficiently.  The function
#' syntax and behavior is largely modeled after that of the
#' \code{cosine()} function from the \code{lsa} package, although
#' with a very different implementation.
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
#' @param inverse
#' Logical; should the inverse covariance matrix be returned?
#' 
#' @return
#' The \eqn{n\times n} matrix of all pair-wise vector cosine
#' similarities of the columns.
#' 
#' @examples
#' x <- matrix(rnorm(10*3), 10, 3)
#' 
#' coop::cosine(x)
#' coop::cosine(x[, 1], x[, 2])
#' 
#' @author Drew Schmidt
#' @seealso \code{\link{sparsity}}
#' @name cosine
#' @rdname cosine
NULL

#' @rdname cosine
#' @export
cosine <- function(x, y, use="everything", inverse=FALSE) UseMethod("cosine")



#' @export
cosine.matrix <- function(x, y, use="everything", inverse=FALSE)
{
  co_matrix(x, y, CO_SIM, use, inplace=FALSE, inverse=inverse)
}

#' @export
cosine.data.frame <- function(x, y, use="everything", inverse=FALSE)
{
  x_mat = as.matrix(x)
  if (!missing(y) && is.data.frame(y))
    y = as.matrix(y)
  
  cosine.matrix(x_mat, y, use=use, inverse=inverse)
}



#' @export
cosine.default <- function(x, y, use="everything", inverse=FALSE)
{
  co_vecvec(x, y, CO_SIM, use)
}



#' @export
cosine.simple_triplet_matrix <- function(x, y, use="everything", inverse=FALSE)
{
  if (!missing(y))
    stop("argument 'y' can not be used with a matrix 'x'")
  
  n <- x$ncol
  a <- x$v
  i <- x$i
  j <- x$j
  index <- 1L
  type <- CO_SIM
  
  if (length(a) != length(i) || length(i) != length(j))
    stop("Malformed simple_triplet_matrix: lengths of 'v', 'i', and 'j' do not agree")
  
  ret <- co_sparse(n, a, i, j, index, type, use, inverse)
  if (!is.null(colnames(x)))
  {
    rownames(ret) <- colnames(x)
    colnames(ret) <- colnames(x)
  }
  
  ret
}



#' @export
cosine.dgCMatrix <- function(x, y, use="everything", inverse=FALSE)
{
  if (!missing(y))
    stop("argument 'y' can not be used with a matrix 'x'")
  
  n <- ncol(x)
  a <- x@x
  i <- x@i
  j <- csc_to_coo(i, x@p)
  index <- 0L
  type <- CO_SIM
  
  if (length(a) != length(i) || length(i) != length(j))
    stop("Malformed dgCMatrix: lengths of 'x', 'i', and 'p' do not agree")
  
  ret <- co_sparse(n, a, i, j, index, type, use, inverse)
  if (!is.null(colnames(x)))
  {
    rownames(ret) <- colnames(x)
    colnames(ret) <- colnames(x)
  }
  
  ret
}



### tcosine

#' @rdname cosine
#' @export
tcosine <- function(x, y, use="everything", inverse=FALSE) UseMethod("tcosine")



#' @export
tcosine.matrix <- function(x, y, use="everything", inverse=FALSE)
{
  co_matrix(x, y, CO_SIM, use, inplace=FALSE, trans=TRUE, inverse=inverse)
}

#' @export
tcosine.data.frame <- function(x, y, use="everything", inverse=FALSE)
{
  x_mat = as.matrix(x)
  if (!missing(y) && is.data.frame(y))
    y = as.matrix(y)
  
  tcosine.matrix(x_mat, y, use=use, inverse=inverse)
}
