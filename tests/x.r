library(coop)
set.seed(1234)

m <- 30
n <- 5
x <- matrix(rnorm(m*n), m, n)

t1 <- cor(x)
t2 <- pcor(x, inplace=TRUE)

t1
t2
all.equal(t1, t2)
