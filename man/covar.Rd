% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covar.r
\name{covar}
\alias{covar}
\alias{tcovar}
\title{Covariance}
\usage{
covar(x, y, use = "everything", inplace = FALSE, inverse = FALSE)

tcovar(x, y, use = "everything", inplace = FALSE, inverse = FALSE)
}
\arguments{
\item{x}{A numeric dataframe/matrix or vector.}

\item{y}{A vector (when \code{x} is a vector) or missing (blank) when 
\code{x} is a matrix.}

\item{use}{The NA handler, as in R's \code{cov()} and \code{cor()}
functions.  Options are "everything", "all.obs", and 
"complete.obs".}

\item{inplace}{Logical; if \code{TRUE} then the method used is slower but
uses less memory than if \code{FALSE}.  See \code{?coop-package}
for details.}

\item{inverse}{Logical; should the inverse covariance matrix be returned?}
}
\value{
The covariance matrix.
}
\description{
An optimized, efficient implemntation for computing covariance.
}
\details{
See \code{?coop-package} for implementation details.
}
\examples{
x <- matrix(rnorm(10*3), 10, 3)

coop::pcor(x)
coop::pcor(x[, 1], x[, 2])

}
\seealso{
\code{\link{cosine}}
}
\author{
Drew Schmidt
}
