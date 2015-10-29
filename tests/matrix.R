if (require(lsa))
{
  x <- matrix(rnorm(30), 10)
  
  t1 <- lsa::cosine(x)
  t2 <- fastcosim::cosine(x)
  stopifnot(all.equal(t1, t2, check.attributes=FALSE))
  
  
  
  x <- t(x)
  
  t1 <- lsa::cosine(x)
  t2 <- fastcosim::cosine(x)
  stopifnot(all.equal(t1, t2, check.attributes=FALSE))
  
  
  
  x <- tcrossprod(x)
  
  t1 <- lsa::cosine(x)
  t2 <- fastcosim::cosine(x)
  stopifnot(all.equal(t1, t2, check.attributes=FALSE))
}


