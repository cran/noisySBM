% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VEMalgorithm.R
\name{mainVEM_Q_par}
\alias{mainVEM_Q_par}
\title{main function of VEM algorithm for fixed number of latent blocks in parallel computing}
\usage{
mainVEM_Q_par(s, ListOfTauRho, modelFamily, model, data, directed)
}
\arguments{
\item{s}{indice of initial point in ListOfTauRho to be used for this run}

\item{ListOfTauRho}{a list of initial points}

\item{modelFamily}{probability distribution for the edges. Possible values:
\code{Gauss}, \code{Gamma}}

\item{model}{Implemented models:
\describe{
  \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
  \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
  \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
  \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
  \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
  \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
  \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
  \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
  \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
  \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
}}

\item{data}{data vector in the undirected model, data matrix in the directed model}

\item{directed}{booelan to indicate whether the model is directed or undirected}
}
\value{
list of estimated model parameters and a node clustering; like the output of fitNSBM()
}
\description{
runs the VEM algorithm the provided initial point
}
