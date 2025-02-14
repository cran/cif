% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mls.R
\name{mls}
\alias{mls}
\title{Multivariate Least-Squares regression}
\usage{
mls(y, x, df_flag = FALSE)
}
\arguments{
\item{y}{left hand side data matrix (one or more columns)}

\item{x}{right hand side data matrix (one or more columns)}

\item{df_flag}{flag = TRUE for degrees of freedom correction for the variance}
}
\value{
out regression coefficients and related statistics
}
\description{
Multivariate Least-Squares regression y = x beta + u
}
\references{
Berta et al. 2020
}
\author{
P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
}
