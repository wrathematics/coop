x <- rnorm(30)
y <- rnorm(30)

t1 <- as.vector(lsa::cosine(x, y))
t2 <- fastcosim::cosine(x, y)
stopifnot(all.equal(t1, t2, check.attributes=FALSE))

