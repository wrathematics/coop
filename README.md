# fastco

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/fastco.png)](https://travis-ci.org/wrathematics/fastco)
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


A for computing covariance, correlation, cosine similarity very quickly.  If you can do this faster, I'd love to know how.

The package is optimized for high performance, and has different implementations for dense matrix inputs, dense vector-vector inputs, and sparse matrix inputs.

For more information, including algorithmic details and benchmarks, see the package vignette.



## Installation

To install the R package:

```r
devtools::install_github("wrathematics/fastco")
```

The source code is also separated from the necessary R wrapper
code.  So it easily builds as a shared library after removing
`src/wrapper.c`.
