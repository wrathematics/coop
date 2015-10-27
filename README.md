# fastcosim

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/fastcosim.png)](https://travis-ci.org/wrathematics/fastcosim)
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


A micro-package for computing cosine similarity of a dense matrix
quickly.  If you can do this faster, I'd love to know how.

The performance should scale fairly well, and will use multiple
threads (if your compiler supports OpenMP) when the matrix has 
more than 2500 columns.  You will see the biggest performance
improvements, in decreasing order of value, by using:

1. Good BLAS.
2. A compiler supporting OpenMP (preferably version 4 or better).

If you're on Linux, you can very easily use OpenBLAS with R.  Users
on other platforms might consider using Revolution R Open, which
ships with Intel MKL.



## The Algorithm

Given an `m`x`n` matrix `x` (input) and an `n`x`n` matrix `cos`
(preallocated output):

1. Compute the upper triangle of the crossproduct `cos = t(x) %*% x` using a symmetric rank-k update (the `_syrk` BLAS function).
2. Iterate over the upper triangle of `cos`:
    1. Divide its off-diagonal values by the square root of the product of its `i`'th and `j`'th diagonal entries.
    2. Replace its diagonal values with 1.
3. Copy the upper triangle of `cos` onto its lower triangle.

The total number of floating point operations is:

1. `m*n*(n+1)` for the symmetric rank-k update.
2. `3/2*(n+1)*n` for the rescaling operation.

The algorithm is `O(mn^2)`, and is dominated by the symmetric rank-k update.



## Benchmarks

Compared to the version in the lsa package (as of 27-Oct-2015),
this implementation performs quite well:

```r
library(rbenchmark)
reps <- 100
cols <- c("test", "replications", "elapsed", "relative")

m <- 2000
n <- 200
x <- matrix(rnorm(m*n), m, n)

benchmark(fastcosim::cosine(x), lsa::cosine(x), columns=cols, replications=reps)

##                   test replications elapsed relative
## 1 fastcosim::cosine(x)          100   0.177    1.000
## 2       lsa::cosine(x)          100 113.543  641.486
```

All benchmarks were performed using:

* R 3.2.2
* OpenBLAS
* gcc 5.2.1
* 4 cores of a Core i5-2500K CPU @ 3.30GHz



## Installation

```r
devtools::install_github("wrathematics/fastcosim")
```

