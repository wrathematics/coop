library(rbenchmark)
reps <- 100
cols <- c("test", "replications", "elapsed", "relative")

m <- 2000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(fastco::cosine(x), lsa::cosine(x), columns=cols, replications=reps)
