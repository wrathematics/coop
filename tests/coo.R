if(require(slam))
{
  library(slam)
  library(fastcosim)
  set.seed(1234)
  
  generate <- function(m, n, size)
  {
    x <- matrix(0, m, n)
    x[sample(m*n, size=size)] <- 10
    x
  }
  
  
  ### Very sparse, has column of 0's
  x <- generate(30, 10, 25)
  coo <- as.simple_triplet_matrix(x)
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
  
  #t1[which(is.nan(t1))] <- 0
  #t2[which(is.nan(t2))] <- 0
  #stopifnot(all.equal(t1, t2))
  
  ### Not very sparse
  x <- generate(30, 10, 75)
  coo <- as.simple_triplet_matrix(x)
  
  t1 <- cosine(x)
  t2 <- cosine(coo)
  
  stopifnot(all.equal(t1, t2))
  
}

