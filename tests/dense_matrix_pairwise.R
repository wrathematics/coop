library(coop)
set.seed(1234)

n <- 21
x = matrix(rnorm(n), 3)
x[sample(n, size=6)] = NA 

t1 <- cov(x, use="pair")
t2 <- covar(x, use="pair")

stopifnot(all.equal(t1, t2))

t1 <- cor(x, use="pair")
t2 <- pcor(x, use="pair")

stopifnot(all.equal(t1, t2))
