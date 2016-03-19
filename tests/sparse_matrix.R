if (require(slam))
{
  library(slam)
  library(coop)
  set.seed(1234)
  
  m <- 30
  n <- 10
  
  ### Very sparse, has column of 0's
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.05)
  coo <- as.simple_triplet_matrix(x)
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
  
  #t1[which(is.nan(t1))] <- 0
  #t2[which(is.nan(t2))] <- 0
  #stopifnot(all.equal(t1, t2))
  
  ### Not very sparse
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.25)
  coo <- as.simple_triplet_matrix(x)
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
}



if (require(Matrix))
{
  library(Matrix)
  library(coop)
  set.seed(1234)
  
  m <- 30
  n <- 10
  
  ### Very sparse, has column of 0's
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.05)
  coo <- as(x, "sparseMatrix")
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
  
  ### Not very sparse
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.25)
  coo <- as(x, "sparseMatrix")
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
}
