library(fastcosim)
library(rbenchmark)

m <- 2000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(cor(x), pcor(x), replications=100)
