library(rbenchmark)
reps <- 100
cols <- c("test", "replications", "elapsed", "relative")

n <- 1000000
x <- rnorm(n)
y <- rnorm(n)

benchmark(fastco::cosine(x, y), lsa::cosine(x, y), columns=cols, replications=reps)
