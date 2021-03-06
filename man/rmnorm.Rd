% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_functions_ramses.R
\name{rmnorm}
\alias{rmnorm}
\title{Random Multivariate Normal Generator}
\usage{
rmnorm(mu, cov, prec)
}
\arguments{
\item{mu}{mean vector.}

\item{cov}{covariance matrix.}

\item{prec}{precision matrix.}
}
\value{
A randomly generated vector with the same length as mu.
}
\description{
Generates a normally distributed vector given a mean vector and
either a covariance matrix or a precision matrix. Does not support
generating multiple vectors at once.
}
