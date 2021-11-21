check_use <- function(use)
{
  match.arg(
    tolower(use),
    c("everything", "all.obs", "complete.obs", "pairwise.complete.obs")
  )
}



is.vec <- function(x)
{
  is.vector(x) && !is.list(x)
}



check.is.flag <- function(x)
{
  if (!(is.logical(x) && length(x) == 1 && (!is.na(x))))
  {
    nm <- deparse(substitute(x))
    stop(paste0("argument '", nm, "' must be TRUE or FALSE"))
  }
  
  invisible()
}



#' @useDynLib coop R_check_badvals
check_badvals <- function(x)
{
  .Call(R_check_badvals, x)
}
