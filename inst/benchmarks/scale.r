library(coop)
library(rbenchmark)
cols <- cols <- c("test", "replications", "elapsed", "relative")
reps <- 25

m <- 10000
n <- 250
x <- matrix(rnorm(m*n, mean=10, sd=3), m, n)

benchmark(scale(x), scaler(x), replications=reps, columns=cols)
