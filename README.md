# fastcosim

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/fastcosim.png)](https://travis-ci.org/wrathematics/fastcosim)
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


A micro-package for computing cosine similarity of a matrix quickly.

Only need this for a homework, so don't expect much more than you see.



## Benchmarks

Compare to the version in the lsa package (as of 26-Oct-2015):

```r
library(rbenchmark)
reps <- 100
cols <- c("test", "replications", "elapsed", "relative")

m <- 2000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(fastcosim::cosine(x), lsa::cosine(x), columns=cols, replications=reps)

##                   test replications elapsed relative
## 1 fastcosim::cosine(x)          100   0.178    1.000
## 2       lsa::cosine(x)          100 113.268  636.337
```

* R 3.2.2
* OpenBLAS
* gcc 5.2.1
* 4 cores of a Core i5-2500K CPU @ 3.30GHz


## Installation

```r
devtools::install_github("wrathematics/fastcosim")
```

