#' @useDynLib coop R_co_mat
#' @useDynLib coop R_co_matmat
#' @useDynLib coop R_co_mat_pairwise
co_matrix <- function(x, y, type, use, inplace, trans=FALSE, inverse=FALSE)
{
  check.is.flag(inplace)
  check.is.flag(inverse)
  if (type != CO_SIM && (inplace && trans))
    stop("Not yet implemented for inplace=TRUE, trans=TRUE, method != cosine()")
  
  
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!missing(y))
    stop("argument 'y' can not be used with a matrix 'x'")
  
  use <- check_use(use)
  if (use == "everything")
  {}
  else if (use == "all.obs")
  {
    if (anyNA(x))
      stop("missing observations in covar/pcor/cosine")
  }
  else if (use == "complete.obs")
  {
    if (anyNA(x))
      x <- naomit(x)
  }
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (use == "pairwise.complete.obs")
    ret <- .Call(R_co_mat_pairwise, x, as.integer(type), as.integer(inverse))
  else
    ret <- .Call(R_co_mat, x, as.integer(type), as.integer(inplace), as.integer(trans), as.integer(inverse))
  
  if (!isTRUE(trans))
  {
    if (!is.null(colnames(x)))
    {
      rownames(ret) <- colnames(x)
      colnames(ret) <- colnames(x)
    }
  }
  else
  {
    if (!is.null(rownames(x)))
    {
      rownames(ret) <- rownames(x)
      colnames(ret) <- rownames(x)
    }
  }
  
  ret
}



#' @useDynLib coop R_naomit_vecvec
#' @useDynLib coop R_co_vecvec
co_vecvec <- function(x, y, type, use)
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
  
  # unlike the matrix version, this should come after casting
  # because I don't feel like doing all this garbage for ints
  use <- match.arg(tolower(use), c("everything", "all.obs", "complete.obs"))
  if (use == "everything")
  {}
  else if (use == "all.obs")
  {
    if (anyNA(x) || anyNA(y))
      stop("missing observations in covar/pcor/cosine")
  }
  else if (use == "complete.obs")
  {
    # perhaps a little hacky...
    out <- .Call(R_naomit_vecvec, x, y)
    x <- out[[1]]
    dim(x) <- NULL
    y <- out[[2]]
    dim(y) <- NULL
  }
  
  .Call(R_co_vecvec, x, y, as.integer(type))
}
