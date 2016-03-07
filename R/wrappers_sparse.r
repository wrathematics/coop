co_sparse <- function(x, y, type)
{
  if (!missing(y))
    stop("argument 'y' can not be used with a matrix 'x'")
  
  a <- x$v
  i <- x$i
  j <- x$j
  n <- as.integer(x$ncol)
  
  if (length(a) != length(i) || length(i) != length(j))
    stop("Malformed simple_triplet_matrix: lengths of 'v', 'i', and 'j' do not agree")
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  if (!is.integer(i))
    storage.mode(i) <- "integer"
  if (!is.integer(j))
    storage.mode(j) <- "integer"
  
  .Call(R_co_sparse, n, a, i, j, as.integer(type))
}
