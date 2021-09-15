library(coop)

f = function(x, y, trans, inv)
{
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"

  .Call("R_co_matmat", x, y, coop:::CO_VAR, FALSE, trans, inv)
}

x = matrix(rnorm(30), 10)
y = matrix(rnorm(30), 10)

a = f(x, y, F, F)
b = cov(x,y)
all.equal(a, b)

a = f(x, y, T, F)
b = cov(t(x), t(y))
all.equal(a, b)



# m = 1000
# n = 250
# 
# x = matrix(rnorm(m*n), m, n)
# y = matrix(rnorm(m*n), m, n)
# 
# library(rbenchmark)
# 
# benchmark( f(x, y, T, F) )






g = function(x, y, trans, inv)
{
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"

  .Call("R_co_matmat", x, y, coop:::CO_ORR, FALSE, trans, inv)
}

x = matrix(rnorm(30), 10)
y = matrix(rnorm(30), 10)

a = g(x, y, F, F)
b = cor(x,y)
all.equal(a, b)

a = g(x, y, T, F)
b = cor(t(x), t(y))
all.equal(a, b)





h = function(x, y, trans, inv)
{
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"

  .Call("R_co_matmat", x, y, coop:::CO_SIM, FALSE, trans, inv)
}

set.seed(1234)
x = matrix(rnorm(30), 10)
y = matrix(rnorm(30), 10)

a = h(x, x, F, F)
b = cosine(x)
all.equal(a, b)

# a = h(x, y, T, F)
# b = cosine(t(x), t(y))
# all.equal(a, b)
