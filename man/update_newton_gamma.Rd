% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VEMalgorithm.R
\name{update_newton_gamma}
\alias{update_newton_gamma}
\title{Perform one iteration of the Newton-Raphson to compute the MLE of the parameters of the Gamma distribution}
\usage{
update_newton_gamma(param, L, M)
}
\arguments{
\item{param}{current parameters of the Gamma distribution}

\item{L}{weighted mean of log(data)}

\item{M}{weighted mean of the data}
}
\value{
updated parameters of the Gamma distribution
}
\description{
Perform one iteration of the Newton-Raphson to compute the MLE of the parameters of the Gamma distribution
}
