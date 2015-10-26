# fastcosim

* **Version:** 0.1-0
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt

A micro-package for computing cosine similarity of a matrix quickly.

Only need this for a homework, so don't expect much more than you see.

## Benchmarks

Compare to the version in the lsa package (as of 26-Oct-2015):

```r
library(rbenchmark)
cols <- c("test", "replications", "elapsed", "relative")

m <- 1000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(fastcosim::cosine(x), lsa::cosine(x), columns=cols)

##                   test replications elapsed relative
## 1 fastcosim::cosine(x)          100   0.404    1.000
## 2       lsa::cosine(x)          100  62.877  155.636
```
