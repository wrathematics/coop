if(require(slam))
{
  library(slam)
  library(fastco)
  set.seed(1234)
  
  m <- 30
  n <- 10
  
  ### Very sparse, has column of 0's
  x <- fastco:::dense_stored_sparse_mat(m, n, prop=.05)
  coo <- as.simple_triplet_matrix(x)
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
  
  #t1[which(is.nan(t1))] <- 0
  #t2[which(is.nan(t2))] <- 0
  #stopifnot(all.equal(t1, t2))
  
  ### Not very sparse
  x <- fastco:::dense_stored_sparse_mat(m, n, prop=.25)
  coo <- as.simple_triplet_matrix(x)
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
  
}
