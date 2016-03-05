cosine_R <- function(x)
{
  cp <- crossprod(x)
  rtdg <- sqrt(diag(cp))
  cos <- cp / tcrossprod(rtdg)
  return(cos)
}

library(compiler)
cosine_R <- cmpfun(cosine_R)

library(fastco)
library(rbenchmark)

m <- 2000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(cosine_R(x), cosine(x), replications=100)
