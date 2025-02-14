% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ec.datadet2.R
\name{ec.datadet2}
\alias{ec.datadet2}
\title{prepares deterministics D^(2)}
\usage{
ec.datadet2(det1, booseas = NA, pntdates = NA)
}
\arguments{
\item{det1}{is the det term with constant and trend created by ec.datadet1.R}

\item{booseas}{is a boolean for daily seasonal dummies}

\item{pntdates}{is a vector of integers where the point dummies should be}
}
\value{
det2mat a matrix with the following columns (daily_seas, point_dummies) and n+npred rows
}
\description{
Prepares deterministic dummies for de-meaned daily seasonal and difference point dummies
}
\references{
Berta et al. 2020
}
\author{
P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
}
