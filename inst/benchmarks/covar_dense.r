library(coop)
library(rbenchmark)
cols <- cols <- c("test", "replications", "elapsed", "relative")
reps <- 25

m <- 10000
n <- 250
x <- matrix(rnorm(m*n), m, n)

cat("### regular\n")
benchmark(cov(x), covar(x), replications=reps, columns=cols)

cat("### pairwise\n")
benchmark(cov(x, use='pair'), covar(x, use='pair'), replications=reps, columns=cols)

