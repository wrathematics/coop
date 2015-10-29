# fastcosim

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/fastcosim.png)](https://travis-ci.org/wrathematics/fastcosim)
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


A micro-package for computing cosine similarity quickly.
If you can do this faster, I'd love to know how.

The package has separate routines for dense matrices/vectors and
for a sparse matrix (like a term-document/document-term matrix).
The use of each is seamless to the user by way of R's S3 methods.

For dense implementations, 
the performance should scale fairly well, and will use multiple
threads (if your compiler supports OpenMP) when the matrix has 
more than 2500 columns.  You will see the biggest performance
improvements, in decreasing order of value, by using:

1. Good BLAS.
2. A compiler supporting OpenMP (preferably version 4 or better).

If you're on Linux, you can very easily use OpenBLAS with R.  Users
on other platforms might consider using Revolution R Open, which
ships with Intel MKL.



## The Algorithms with Notes on Implementation

#### Dense Matrix Input

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

The algorithmic complexity is `O(mn^2)`, and is dominated by the symmetric rank-k update.

#### Dense Vector-Vector Input

Given two `n`-length vectors `x` and `y` (inputs):

1. Compute `crossprod = t(x) %*% y` (using the `_gemm` BLAS function).
2. Compute the square of the Euclidean norms of `x` and `y` (using the `_syrk` BLAS function).
3. Divide `crossprod` from 1 by the square root of the product of the norms from 2.

The total number of floating point operations is:

1. `2n-1` for the crossproduct.
2. `4*n-2` for the two (square) norms.
3. `3` for the division and square root/product.

The algorithmic complexity is `O(n)`.


#### Sparse Matrix Input

Given a matrix stored as a COO with row/column indices `i` and `j`
**where they are sorted by columns first, then rows**, and
corresponding data `a` and number of columns `n` (inputs), and a
preallocated `n`x`n` dense matrix `cos` (output):

1. Initialize `cos` to 0.
2. For each column `j` of `a` (call it `x`), find its first and final position in the COO storage.
    i. If `x` is missing (its entries are all 0), set the `j`'th row and column of the lower triangle of `cos` to 0.  Go to 2.
    ii. Otherwise, for each column `i>j` of `a` (call it `y`), find its first and final position  in the COO storage.
    iii. Compute the dot product of `x` and `y`, `xy`.
    iv. If the dot product is greater than epsilon (`1e-10` for us):
        - Compute the dot products of `x` with itself `xx` and `y` with itself `yy`.
        - Set the `(i, j)`'th entry of `cos` to `xy`/`sqrt(xx*yy)`.
3. Copy the lower triangle to the upper and set the diagonal to 1.



## Benchmarks

All benchmarks were performed using:

* R 3.2.2
* OpenBLAS
* gcc 5.2.1
* 4 cores of a Core i5-2500K CPU @ 3.30GHz
* Linux kernel 4.2.0-16

#### Dense Matrix Input

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

#### Dense Vector-Vector Input

Here the two perform identically:

```r
library(rbenchmark)
reps <- 100
cols <- c("test", "replications", "elapsed", "relative")

n <- 1000000
x <- rnorm(n)
y <- rnorm(n)

benchmark(fastcosim::cosine(x, y), lsa::cosine(x, y), columns=cols, replications=reps)

##                      test replications elapsed relative
## 1 fastcosim::cosine(x, y)          100   0.757    1.000
## 2       lsa::cosine(x, y)          100   0.768    1.015
```


#### Sparse Matrix Input

TODO



## Installation

```r
devtools::install_github("wrathematics/fastcosim")
```

