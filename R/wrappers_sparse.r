co_sparse <- function(n, a, i, j, index, type)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  if (!is.integer(i))
    storage.mode(i) <- "integer"
  if (!is.integer(j))
    storage.mode(j) <- "integer"
  
  .Call(R_co_sparse, as.integer(n), a, i, j, as.integer(index), as.integer(type))
}



csc_to_coo <- function(row_ind, col_ptr)
{
  .Call(R_csc_to_coo, row_ind, col_ptr)
}
