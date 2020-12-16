    #' simulation of a graph according the noisy stochastic block model
    #'
    #' @param n number of nodes
    #' @param theta model parameters of the noisy stochastic block model
    #' \describe{
    #'   \item{pi}{latent block proportions, Q-vector}
    #'   \item{w}{connectivity parameters, N_Q-vector}
    #'   \item{nu0}{parameters of the null distribution}
    #'   \item{nu}{parameters of the alternative distribution}
    #' }
    #' @param modelFamily probability distribution for the edges. Possible values:
    #'       \code{Gauss}, \code{Gamma}, \code{Poisson}
    #' @param directed  indicates if the graph is directed (boolean)
    #'
  #' @return a list with:
  #' \describe{
  #'   \item{dataMatrix}{simulated matrix from the noisy stochastic block model}
  #'   \item{theta}{model parameters of the noisy stochastic block model}
  #'   \item{latentZ}{underlying latent node memberships}
  #'   \item{latentAdj}{underlying latent binary graph}
  #' }
  #' @export
  #'
  #' @examples
  #' n <- 10
  #' Q <- 2
  #' theta <- list(pi= rep(1/Q,Q), nu0=c(0,1))
  #' theta$nu <- matrix(c(-2,10,-2, 1,1,1),nrow=Q*(Q+1)/2,ncol=2)
  #' theta$w <- c(.5, .9, .3)
  #' obs <- rnsbm(n, theta, modelFamily='Gauss')
  #' obs
  rnsbm <-  function(n, theta, modelFamily='Gauss', directed=FALSE){
    N <- if (directed) n*(n-1) else n*(n-1)/2
    Q <- length(theta$pi)

    # latent variables
    Z <- sample(1:Q, n, replace=TRUE, prob=theta$pi)

    # adjacency matrix
    A <- matrix(0, n, n)
    if (directed){
      for (i in 1:n){
        A[i,-i] <- stats::rbinom(n-1, 1, theta$w[convertGroupPair(Z[i], Z[-i], Q, directed)])
      }
    }else{
      for (i in 1:(n-1)){
        A[i,(i+1):n] <- stats::rbinom(n-i, 1, theta$w[convertGroupPair(Z[i],Z[(i+1):n], Q, directed)])
      }
    }

    # noisy observations under the null
    X <- matrix(0, n, n)
    nonZero <- if (directed) !(diag(n)) else upper.tri(diag(n))
    if (modelFamily=='Gauss'){
      X[nonZero] <- stats::rnorm(N, theta$nu0[1], theta$nu0[2])
    }
    if (modelFamily=='Gamma'){
      X[nonZero] <- stats::rgamma(N, theta$nu0[1], theta$nu0[2])
    }
    if (modelFamily=='Poisson'){
      X[nonZero] <- stats::rpois(N, theta$nu0)
    }

    m <- if (directed) n else n-1
    for (i in 1:m){
      nonzeroind <- which(A[i,]!=0)
      L <- length(nonzeroind)
      if (L>=1){
        if (modelFamily=='Gauss'){
          ind_i <- convertGroupPair(Z[i], Z[nonzeroind], Q, directed)
          X[i, nonzeroind] <- stats::rnorm(L, theta$nu[ind_i,1], theta$nu[ind_i,2])
        }
        if (modelFamily=='Gamma'){
          ind_i <- convertGroupPair(Z[i], Z[nonzeroind], Q, directed)
          X[i, nonzeroind] <- stats::rgamma(L, theta$nu[ind_i,1], theta$nu[ind_i,2])
        }
        if (modelFamily=='Poisson'){
          X[i, nonzeroind] <- stats::rpois(L, theta$nu[convertGroupPair(Z[i],Z[nonzeroind],Q,directed)])
        }
      }
    }

    if (!directed){
      A <- A + t(A)
      X <- X + t(X)
    }
    return(list(dataMatrix=X, theta=theta, latentZ=Z, latentAdj=A))
  }

