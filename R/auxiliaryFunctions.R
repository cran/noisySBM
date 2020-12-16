#' transform a pair of nodes (i,j) into an identifying integer
#'
#' Associates an identifying integer with a pair of nodes (i,j)
#'
#' @details returns the row number of the matrix build by listNodePairs(n)
#'     containing the pair (i,j)
#'
#' @param i scalar or vector
#' @param j scalar or vector, same length as i
#' @param n number of vertices
#' @param directed booelan to indicate whether the model is directed or undirected
convertNodePair <- function(i,j,n, directed){
  if (sum((i>n) | (j>n))>0){
    stop("Your index is out of range")
  }
  if (directed){#directed case
    dyads <- (i-1)*(n-1)+j-(i<j)
  } else {#undirected case
    dyads <- c(0,cumsum((n-1):1))[pmin(i,j)] + abs(j-i)
  }
  return(dyads)
}



#' returns a list of all possible node pairs (i,j)
#'
#' @param n number of nodes
#' @param directed indicates if the graph is directed
#'
#' @return a 2-column matrix with all possible node pairs (i,j)
listNodePairs <- function(n, directed=FALSE){
  if (!exists('N'))
    N <- if (directed) n*(n-1) else n*(n-1)/2
  index <- matrix(0,N,2)
  if (directed){ # directed
    index[,1] <- rep(1:n,each=n-1)
    k <- (1:n^2)[-seq(1,n^2,by=n+1)]
    index[,2] <- rep(1:n,n)[k]
  }else { # undirected
    index[,1] <- rep(1:(n-1),times=(n-1):1)
    toto <- c()
    for (k in 1:(n-2)){
      toto <- c(toto, k*(n-1)+ 1:k)
    }
    index[,2] <- rep(2:n,n-1)[-toto]
  }
  return(index)
}




#' transform a pair of block identifiers (q,l) into an identifying integer
#'
#' this is the inverse function of convertGroupPairIdentifier()
#'
#' @param q indicator of a latent block
#' @param l indicator of a latent block
#' @param Q number of latent blocks
#' @param directed indicates if the graph is directed
convertGroupPair <- function(q, l, Q, directed){
  if (directed){
    index <- (q-1)*Q+l
  } else { # undirected
    qp <- pmin(q,l)
    lp <- pmax(q,l)
    index <- (2*Q-qp+2)*(qp-1)/2 +lp-qp+1
  }
  return(index)
}


#' takes a scalar indice of a group pair (q,l) and returns the values q and l
#'
#' this is the inverse function of convertGroupPair()
#'
#' @param ind_ql indicator for a pair of latent blocks
#' @param Q number of latent blocks
convertGroupPairIdentifier <- function(ind_ql, Q){
  w <- cumsum((Q-1):1)
  q <- which.max(ind_ql<=w)
  w <- c(0, w)
  l <- ind_ql - w[q] + q
  return(c(q,l))
}


#' corrects values of the variational parameters tau that are too close to the 0 or 1
#'
#' @param tau variational parameters
correctTau <- function(tau){
  tau <- pmin(tau,.Machine$double.xmax)
  tau <- pmax(tau,.Machine$double.xmin)
  tau <- tau/sum(tau)
  tau <- pmin(tau,1-1e-7)
  tau <- pmax(tau,1e-7)
  tau <- tau/sum(tau)

  return(tau)
}


#' evaluate the density in the current model
#'
#' @param x vector with points where to evaluate the density
#' @param nu distribution parameter
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}, \code{Poisson}
modelDensity <- function(x, nu, modelFamily='Gauss'){
  if (modelFamily=='Gauss')
    res <- stats::dnorm(x, nu[1], nu[2])
  if (modelFamily=='Gamma')
    res <- stats::dgamma(x, nu[1], nu[2])
  res[res<=.Machine$double.eps] <- .Machine$double.eps
  return(res)
}


#' Evaluate tau_q*tau_l in the noisy stochastic block model
#'
#' @param q indicator of a latent block
#' @param l indicator of a latent block
#' @param tau variational parameters
#' @param n number of vertices
#' @param directed booelan to indicate whether the model is directed or undirected
getTauql <- function(q, l, tau, n, directed){
  n <- ncol(tau)
  Q <- nrow(tau)
  N <- if (directed) n*(n-1) else n*(n-1)/2

  ind.all <- listNodePairs(n, directed)
  if (Q==1)
    tauql <- rep(1, N)
  else{
    if ((q==l))
      tauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]]
    else
      tauql <- tau[q, ind.all[,1]]*tau[l, ind.all[,2]] + tau[q, ind.all[,2]]*tau[l, ind.all[,1]]
  }
  return(tauql)
}














