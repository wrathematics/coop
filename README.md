# coop

* **Version:** 0.6-2
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/coop.png)](https://travis-ci.org/wrathematics/coop)
* **License:** [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


The **coop** package does co-operations: covariance, correlation, and cosine similarity.  And it does them very quickly.  If you can do this faster, I'd love to know how.

The package is optimized for high performance, and has different implementations for dense matrix inputs, dense vector-vector inputs, and sparse matrix inputs.  Note that to get good performance with this package (as seen in these benchmarks), you will need to use a good BLAS library.  See the package vignette for details.

For more information, including algorithmic details, see the package vignettes.



## Installation

To install the R package, run:

```r
install.packages("coop")
```

The development version is maintained on GitHub, and can easily be installed by any of the packages that offer installations from GitHub:

```r
### Pick your preference
devtools::install_github("wrathematics/coop")
ghit::install_github("wrathematics/coop")
remotes::install_github("wrathematics/coop")
```

The C internals are completely separated from the R wrapper code.  So if you prefer, you can easily build this as a C shared library after removing the file `src/wrapper.c`.





## Package Use

The package has functions for covariance and pearson correlation with interfaces that mimic base R's, with the addition of a `cosine()` function.  At this time, the basic interface looks like this:

```r
### matrix input
covar(x)     # like cov(x)
pcor(x)      # like cor(x)
cosine(x)    # like lsa::cosine(x)

### vector input
covar(x, y)  # like cov(x, y)
pcor(x, y)   # like cor(x, y)
cosine(x, y) # like lsa::cosine(x, y)
```

There are also `t` versions of the functions which operate on the transposed data (without producing a copy).  So `tcovar(x)` will do the same computation as as `cov(t(x))` (but much more efficiently).

The functions also have an additional argument `inverse`, which will return the matrix inverse of the specified operation.  So `covar(x, inverse=TRUE)` will return the inverted covariance matrix.

For more details, see the package vignette.





## Benchmarks

Here we provide some benchmarks for dense matrices.  The package also has vector-vector methods for each operation, and a sparse method for cosine similarity.  These also perform quite well, but in the case of the former are generally not performance intensive, and in the case of the latter, I am not aware of any other sparse cosine similarity implementations available to R.

All of these benchmarks can be found in the source tree of this package, under `coop/inst/benchmarks`.  Implementation details can be found in the package vignette.  All tests performed using:

* R 3.2.3
* OpenBLAS
* gcc 5.3.1
* 4 cores of a Core i5-2500K CPU @ 3.30GHz

#### Setup

```r
library(coop)
library(rbenchmark)
cols <- c("test", "replications", "elapsed", "relative")
reps <- 25

m <- 10000
n <- 250
x <- matrix(rnorm(m*n), m, n)
```

#### Covariance

```r
benchmark(cov(x), covar(x), replications=reps, columns=cols)
##       test replications elapsed relative
## 2 covar(x)           25   0.431    1.000
## 1   cov(x)           25   8.656   20.084
```

#### Correlation

```r
benchmark(cor(x), pcor(x), replications=reps, columns=cols)
##      test replications elapsed relative
## 1  cor(x)           25   8.695   16.818
## 2 pcor(x)           25   0.517    1.000
```

#### Cosine

```r
benchmark(lsa::cosine(x), cosine(x), replications=reps, columns=cols)
##             test replications elapsed relative
## 2      cosine(x)           25   0.341    1.000
## 1 lsa::cosine(x)           25 181.045  530.924
```

We can outperform the **lsa** implementation using just R code as follows:

```r
cosine_R <- function(x)
{
  cp <- crossprod(x)
  rtdg <- sqrt(diag(cp))
  cos <- cp / tcrossprod(rtdg)
  return(cos)
}
```

We note that while this implementation is reasonably "clock efficient", it is very memory wasteful compared to the implementation in **coop** (and for much larger data sizes, the implementation in **coop** will dominate).

```r
library(compiler)
cosine_R <- cmpfun(cosine_R)

benchmark(cosine_R(x), cosine(x), replications=reps, columns=cols)
##          test replications elapsed relative
## 1 cosine_R(x)           25   0.287    1.079
## 2   cosine(x)           25   0.266    1.000
```

In fact, similar tricks can be played with both covariance and correlation.  One major reason **coop** is so much faster is because of its careful use of the BLAS.  See the vignette for more information.
