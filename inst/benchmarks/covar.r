library(fastco)
library(rbenchmark)

m <- 2000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(cov(x), covar(x), replications=100)
