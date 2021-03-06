% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VEMinitialization.R
\name{getRho}
\alias{getRho}
\title{compute rho associated with given values of w, nu0 and nu}
\usage{
getRho(Q, w, nu0, nu, data, modelFamily)
}
\arguments{
\item{Q}{number of latent blocks in the noisy stochastic block model}

\item{w}{weight parameter in the noisy stochastic block model}

\item{nu0}{null parameter in the noisy stochastic block model}

\item{nu}{alternative parameter in the noisy stochastic block model}

\item{data}{data vector in the undirected model, data matrix in the directed model}

\item{modelFamily}{probability distribution for the edges. Possible values:
\code{Gauss}, \code{Gamma}}
}
\value{
a matrix of conditional probabilities of an edge given the node memberships of the interacting nodes
}
\description{
compute rho associated with given values of w, nu0 and nu
}
