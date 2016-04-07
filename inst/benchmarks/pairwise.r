library(coop)
library(rbenchmark)
cols <- cols <- c("test", "replications", "elapsed", "relative")
reps <- 25

m <- 1000
n <- 150
x <- matrix(rnorm(m*n), m, n)

x[sample(m*n, size=.1*m*n)] <- NA_real_

benchmark(R=(res.r <- cov(x, use="pair")), coop=(res.coop <- covar(x, use="pair")), replications=reps, columns=cols)
all.equal(res.r, res.coop)
