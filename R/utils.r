check_use <- function(use)
{
  match.arg(tolower(use), c("everything", "all.obs", "complete.obs"))
}



is.vec <- function(x)
{
  is.vector(x) && !is.list(x)
}
