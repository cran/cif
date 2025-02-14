% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ec.datadet1.R
\name{ec.datadet1}
\alias{ec.datadet1}
\title{prepares deterministics D^(1)}
\usage{
ec.datadet1(n, befpn, breaks)
}
\arguments{
\item{n}{is the number of obs in available data}

\item{befpn}{is a vector with (begtrim,endtrim,nforecast,npred,nhstar)}

\item{breaks}{is a vector of integers where the trend breaks should be}
}
\value{
matdet1 a matrix with the following columns (1_vec, t_vec) and (n+npred) rows
}
\description{
Prepares deterministic data
}
\references{
Berta et al. 2020
}
\author{
P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
}
