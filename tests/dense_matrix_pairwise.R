library(coop)
set.seed(1234)

x = matrix(rnorm(9), 3)
x[2,1] = x[3,2] = NA 

t1 <- cov(x, use="pair")
t2 <- covar(x, use="pair")

stopifnot(all.equal(t1, t2))
