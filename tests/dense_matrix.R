### Cosine similarity
cosine <- function(x)
{
  cp <- crossprod(x)
  rtdg <- sqrt(diag(cp))
  cos <- cp / tcrossprod(rtdg)
  return(cos)
}

check <- function(x)
{
  t1 <- cosine(x)
  t2 <- coop::cosine(x)
  stopifnot(all.equal(t1, t2, check.attributes=FALSE))
}

x <- matrix(rnorm(30), 10)
check(x)
check(t(x))
check(crossprod(x))
check(tcrossprod(x))



### Pearson correlation
check <- function(x)
{
  t1 <- cor(x)
  t2 <- coop::pcor(x)
  stopifnot(all.equal(t1, t2, check.attributes=FALSE))
}

check(x)
check(t(x))
check(crossprod(x))
check(tcrossprod(x))



### Covariance
check <- function(x)
{
  t1 <- cov(x)
  t2 <- coop::covar(x)
  stopifnot(all.equal(t1, t2, check.attributes=FALSE))
}

check(x)
check(t(x))
check(crossprod(x))
check(tcrossprod(x))
