% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ec.main.R
\name{ec.main}
\alias{ec.main}
\title{Forecast with Vector Error Correction Model}
\usage{
ec.main(
  y,
  ndet = c(2, 1),
  nlag,
  befpn,
  breaks = NA,
  booseas = NA,
  pntdates = NA,
  drop1 = NA,
  drop2 = NA,
  cal,
  kval = 1.959964
)
}
\arguments{
\item{y}{matrix with time across rows and variables in columns}

\item{ndet}{vector of lenght 3, (i,j,q): i for EG1-st stage, j for EG-2nd stage, q number of breaks
i,j=0 no deterministics
i,j=1 constant
i,j=2 constant and trend}

\item{nlag}{number of lags in the VAR}

\item{befpn}{begtrim, endtrim, nforecast, npred}

\item{breaks}{vector with observation numbers for T1,T2,...}

\item{booseas}{boolean =T if seasonal dummies, =F otherwise}

\item{pntdates}{vector with observation numbers for point dummies}

\item{drop1}{selection of det1 regressors in first stage to drop}

\item{drop2}{selection of det1 regressors in second stage to drop}

\item{cal}{calendar for the y matrix}

\item{kval}{how many se to use, default kval=1.959964}
}
\value{
results Output contains the a set of estimates and forecasting results.
}
\description{
This function estimate VECM model. Selects begtrim and entrim period, define lag and run.
}
\references{
Berta et al. 2020
}
\author{
P. Berta, P. Paruolo, S. Verzillo, PG. Lovaglio
}
