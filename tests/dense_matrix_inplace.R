library(coop)
set.seed(1234)

m <- 30
n <- 5
x <- matrix(rnorm(m*n), m, n)

t1 <- cov(x)
t2 <- covar(x, inplace=TRUE)

all.equal(t1, t2)
