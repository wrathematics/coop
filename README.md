# fastcosim

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/fastcosim.png)](https://travis-ci.org/wrathematics/fastcosim)
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


A micro-package for computing cosine similarity very quickly.  If you can do this faster, I'd love to know how.

The package is optimized for high performance, and has different implementations for dense matrix inputs, dense vector-vector inputs, and sparse matrix inputs.

For more information, including algorithmic details and benchmarks, see the package vignette.



## Installation

To install the R package:

```r
devtools::install_github("wrathematics/fastcosim")
```

The source code is also separated from the necessary R wrapper
code.  So it easily builds as a shared library after removing
`src/wrapper.c`.
