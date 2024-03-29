% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal.R
\name{normal_subsample}
\alias{normal_subsample}
\title{Find a subset of points forming multivariate normal distribution.}
\usage{
normal_subsample(matrix, size)
}
\arguments{
\item{matrix}{n by m R matrix where columns are samples and rows are
corresponding components of vectors.}

\item{size}{the number of points to be subsetted.}
}
\description{
The method takes a set of n-dimensional points of size m and tries 
to select \code{size} points from the set by minimizing difference
between emperical characteristic function of subset and characteristic funciton
of multivariate normal with sample mean and covariance from the whole set
as parameters.
}
