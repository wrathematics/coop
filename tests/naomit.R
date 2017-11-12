library(coop)
set.seed(1234)

naomit = coop:::naomit
check = function(a, b) stopifnot(all.equal(a, b, check.attributes=FALSE))

test = function(m, n, prop=.01)
{
  x = matrix(rnorm(m*n, sd=10000), m, n)
  
  check(na.omit(x), naomit(x))
  
  x[sample(m*n, size=m*n*prop)] = NA
  
  truth = na.omit(x)
  if (any(dim(truth) == 0))
    stop("zero rows - bad test seed/prop values")
  
  check(truth, naomit(x))
  
  storage.mode(x) = "integer"
  check(na.omit(x), naomit(x))
}

test(100, 20)
test(10, 2)



### TODO
# if (require(slam))
# {
#   library(slam)
#   csc = as.simple_triplet_matrix(y)
#   z = as.matrix(coop:::naomit_coo(as.double(csc$v), csc$i, csc$j))
#   stopifnot(all.equal(na.omit(y), z, check.attributes=FALSE))
# }
