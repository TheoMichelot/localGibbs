% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/n2w.R
\name{n2w}
\alias{n2w}
\title{Transform parameters from natural to working scale}
\usage{
n2w(beta, rdist = c("fixed", "multistate", "exp", "gamma", "weibull"),
  r = NULL, sigma = NULL, shape = NULL, rate = NULL, gamma = NULL,
  xy = NULL, norm = FALSE)
}
\arguments{
\item{beta}{Coefficients of the RSF}

\item{rdist}{Distribution of the availability radius ("fixed", "multistate", 
"exp", "gamma", or "weibull)}

\item{r}{Availability radius, if rdist="fixed" (or vector of radii if multistate)}

\item{sigma}{Vector of standard deviations (if norm=TRUE)}

\item{shape}{Shape parameter, if rdist="gamma"}

\item{rate}{Rate parameter, if rdist="exp" or rdist="gamma"}

\item{gamma}{Transition probability matrix, if rdist="multistate"}

\item{xy}{Matrix of observed locations, needed to derive maximum step length
if rdist="fixed"}

\item{norm}{Logical. TRUE if normal transition density. (Only for multistate case)}
}
\value{
Vector of parameters on the working scale
}
\description{
Transform parameters from natural to working scale
}
