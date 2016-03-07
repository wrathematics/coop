is.vec <- function(x)
{
  is.vector(x) && !is.list(x)
}



co_matrix <- function(x, y, type)
{
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!missing(y))
    stop("argument 'y' can not be used with a matrix 'x'")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(R_co_mat, x, as.integer(type))
}



co_vecvec <- function(x, y, type)
{
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  
  if (missing(y) && type != CO_VAR)
    return(1.0)
  else if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  else if (!is.vec(y))
    stop("argument 'y' must be a non-list vector")
  
  if (length(x) != length(y))
    stop("vectors 'x' and 'y' must have the same length")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  .Call(R_co_vecvec, x, y, as.integer(type))
}
