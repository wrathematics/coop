if(require(slam))
{
  library(slam)
  library(fastcosim)
  
  set.seed(1234)
  
  m <- 30
  n <- 10
  x <- matrix(0, m, n)
  x[sample(m*n, size=25)] <- 10
  
  coo <- as.simple_triplet_matrix(x)
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  
  stopifnot(all.equal(t1, t2))
  
  #t1[which(is.nan(t1))] <- 0
  #t2[which(is.nan(t2))] <- 0
  #stopifnot(all.equal(t1, t2))
}

