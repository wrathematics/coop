library(coop)

m <- 20
n <- 5
x <- matrix(rnorm(m*n, mean=10, sd=3), m, n)

# center/scale
t1 <- scale(x)
t2 <- scaler(x)
stopifnot(all.equal(t1, t2))

# center
t1 <- scale(x, scale=FALSE)
t2 <- scaler(x, scale=FALSE)
stopifnot(all.equal(t1, t2))

# scale
t1 <- scale(x, center=FALSE)
t2 <- scaler(x, center=FALSE)
stopifnot(all.equal(t1, t2))
