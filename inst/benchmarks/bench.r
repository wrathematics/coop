library(rbenchmark)
reps <- 100
cols <- c("test", "replications", "elapsed", "relative")

m <- 1000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(fastcosim::cosine(x), lsa::cosine(x), columns=cols, replications=reps)
