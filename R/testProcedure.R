#' compute conditional l-values in the noisy stochastic block model
#'
#' @param dataVec data vector
#' @param Z a node clustering
#' @param theta list of parameters for a noisy stochastic block model
#' @param directed indicates if the graph is directed
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss} and \code{Gamma}
#'
#' @return conditional l-values in the noisy stochastic block model
lvaluesNSBM <- function(dataVec, Z, theta, directed=FALSE, modelFamily='Gauss'){
  n <- length(Z)
  Q <- length(theta$pi)

  Z_matrix <- matrix(0, nrow=Q, ncol=n)
  ind <- matrix(c(Z,1:n),ncol=2)
  Z_matrix[ind] <- 1

  lval <- rep(0, length(dataVec))
  ind <- 0
  ind.all <- listNodePairs(n, directed)
  fnu0data <- modelDensity(dataVec, theta$nu0, modelFamily)
  for (q in 1:Q){
    for (l in q:Q){
      ind <- ind + 1
      mask <- if (q==l) (Z_matrix[q, ind.all[,1]]*Z_matrix[l, ind.all[,2]]) else
        (Z_matrix[q, ind.all[,1]]*Z_matrix[l, ind.all[,2]] + Z_matrix[q, ind.all[,2]]*Z_matrix[l, ind.all[,1]] )
      part_f0 <- fnu0data*(1-theta$w[ind])
      lval <- lval + part_f0/ (modelDensity(dataVec, theta$nu[ind,], modelFamily) * theta$w[ind] + part_f0)*mask
    }
  }

  return(lval)
}




#' auxiliary function for the computation of q-values
#'
#' @param theta list of parameters for a noisy stochastic block model
#' @param ind indicator for a pair of latent blocks
#' @param t l-values
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss} and \code{Gamma}
q_delta_ql <- function(theta, ind, t, modelFamily='Gauss'){
  # for a given (q,l)
  if (modelFamily=='Gauss'){
    mu0 <- theta$nu0[1]
    sigma0 <- theta$nu0[2]
    mu <- theta$nu[ind, 1]
    sigma <- theta$nu[ind, 2]

    nb_t <- length(t)
    res <- matrix(0, nb_t, 2)
    res[t==1,] <- 1
    ind_t <- (t>0)&(t<1)
    if (length(ind_t)>0){
      a <- sigma^(-2)-sigma0^(-2) # scalar
      b <- -2*(mu/sigma^2 -mu0/sigma0^2)    # scalar
      cVec <- mu^2/sigma^2 -mu0^2/sigma0^2 + 2*log(sigma/sigma0*(1/theta$w[ind]-1)*(1/t-1))    # vector
      if(a!=0){
        res[ind_t,] <- res[ind_t,] + (a<0)
        ind_2 <- (b^2>(4*a*cVec))&ind_t
        if(sum(ind_2)>0){
          z <- (-b+matrix(sqrt(pmax(0,b^2-4*a*cVec[ind_2])),ncol=1)%*%c(1,-1))/(2*a) # matrix sum(ind_2) x 2
          res[ind_2,2] <- res[ind_2,2] + stats::pnorm(z[,1], mu, sigma) - stats::pnorm(z[,2], mu, sigma)
          res[ind_2,1] <- res[ind_2,1] + stats::pnorm(z[,1], mu0, sigma0) - stats::pnorm(z[,2], mu0, sigma0)
        }
      }else{
        if (b!=0){
          res[ind_t,2] <- if (b<0) 1-stats::pnorm(-cVec[ind_t]/b, mu, sigma) else stats::pnorm(-cVec[ind_t]/b, mu, sigma)
          res[ind_t,1] <- if (b<0) 1-stats::pnorm(-cVec[ind_t]/b, mu0, sigma0) else stats::pnorm(-cVec[ind_t]/b, mu0, sigma0)
        }else{
          res[ind_t,] <- 1*(t[ind_t] >= 1- theta$w[ind])
        }
      }
    }
  }
  if (modelFamily=='Gamma')  {
    nu0 <- theta$nu0
    nu <- theta$nu[ind, ]
    ## warning si nu0[1] ! = nu[1]
    if (nu0[1] != nu[1]){
      stop('shapes between null and alternative distribution are different')
      } else {
      nb_t <- length(t)
      res <- matrix(0, nb_t, 2)
      res[t==1,] <- 1
      ind_t <- (t>0)&(t<1)
      if (length(ind_t)>0){
        if (nu0[2]>nu[2]){
          res[ind_t,1] <- 1- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$w[ind]-1)*(nu0[2]/nu[2])^nu0[1]), nu0[1] , nu0[2])
          res[ind_t,2] <- 1- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$w[ind]-1)*(nu0[2]/nu[2])^nu0[1]), nu[1], nu[2])
        }

        if (nu0[2]< nu[2]){
          res[ind_t,1] <- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$w[ind]-1)*(nu0[2]/nu[2])^nu0[1]),  nu0[1] , nu0[2])
          res[ind_t,2] <- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$w[ind]-1)*(nu0[2]/nu[2])^nu0[1]), nu[1] , nu[2])
        }

        if (nu0[2]==nu[2]){
          res[ind_t,1] <- 1*(1-theta$w[ind] <= t[ind_t])
          res[ind_t,2] <- 1*(1-theta$w[ind] <= t[ind_t])
        }
      }
    }
  }
  return(res)
}




#' compute q-values in the noisy stochastic block model
#'
#' @param dataVec  data vector
#' @param Z a node clustering
#' @param theta list of parameters for a noisy stochastic block model
#' @param lvalues conditional l-values in the noisy stochastic block model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss} and \code{Gamma}
#' @param directed indicates if the graph is directed
#'
#' @return q-values in the noisy stochastic block model
qvaluesNSBM <- function(dataVec, Z, theta, lvalues, modelFamily='Gauss', directed=FALSE){
  Q <- length(theta$pi)
  num <- den <- result <- rep(0, length(dataVec))
  ind <- 0
  ind.all <- listNodePairs(length(Z), directed)
  for (q in 1:Q){
    for (l in q:Q){
      ind <- ind + 1
      q01_lvalues_seq <- theta$pi[q]*theta$pi[l]*q_delta_ql(theta, ind, lvalues, modelFamily)
      f <- 1 + (q!=l)
      num <- num + f*(1-theta$w[ind]) * q01_lvalues_seq[,1]
      den <- den + f*theta$w[ind]*q01_lvalues_seq[,2]
    }
  }
  den <- den + num
  ind <- (den!=0)
  result[ind] <- num[ind]/den[ind]
  return(result)
}

#' new graph inference procedure
#'
#' @details graph inference procedure based on conditional q-values in the noisy stochastic block model. It works in the
#' Gaussian model, and also in the Gamma model, but only if the shape parameters of the
#' Gamma distributions under the null and the alternatives are identical (e.g. when all distributions
#' are exponentials).
#'
#' @param dataMatrix observed adjacency matrix, nxn matrix
#' @param nodeClustering n-vector of hard node Clustering
#' @param theta parameter of the noisy stochastic block model
#' @param alpha confidence level
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss} and \code{Gamma}
#'
#' @return a list with:
#'  \describe{
#'     \item{\code{A}}{resulting binary adjacency matrix}
#'     \item{\code{qvalues}}{vector with conditional q-values in the noisy stochastic block model}
#'  }
#' @export
#'
#' @examples
#' set.seed(1)
#' theta <- list(pi=c(.5,.5), w=c(.8,.1,.2), nu0=c(0,1), nu=matrix(c(-1,5,10, 1,1,1), ncol=2))
#' obs <- rnsbm(n=30, theta)
#' # res_gauss <- fitNSBM(obs$dataMatrix, nbCores=1)
#' resGraph <- graphInference(obs$dataMatrix, res_gauss[[2]]$clustering, theta, alpha=0.05)
#' sum((resGraph$A))/2 # nb of derived edges
#' sum(obs$latentAdj)/2 # correct nb of edges
graphInference <- function(dataMatrix, nodeClustering, theta, alpha=0.05, modelFamily='Gauss'){

  dataVec <- dataMatrix[lower.tri(dataMatrix)]
  directed <- FALSE
  lval_results <- lvaluesNSBM(dataVec, nodeClustering, theta, directed, modelFamily)
  qval_results <- qvaluesNSBM(dataVec, nodeClustering, theta, lval_results, modelFamily, directed)

  n <- length(nodeClustering)
  A <- matrix(0, n, n)
  A[lower.tri(A)] <- (qval_results<alpha)
  A <- A + t(A)
  return(list(A=A, qvalues=qval_results))
}

#' plot the data matrix, the inferred graph and/or the true binary graph
#'
#' @param dataMatrix observed data matrix
#' @param inferredGraph graph inferred by the multiple testing procedure via graphInference()
#' @param binaryTruth true binary graph
#'
#' @return a list of FDR and TDR values, if possible
#' @export
plotGraphs <- function(dataMatrix=NULL, inferredGraph=NULL, binaryTruth=NULL){
  res <- NULL
  if (!is.null(dataMatrix))
    stats::heatmap(dataMatrix, Rowv=NA, Colv=NA, symm=TRUE, col = RColorBrewer::brewer.pal(9,'YlGn'),
            xlab='data matrix')
  xlab2 <- 'inferred graph'
  xlab3 <- 'true binary graph'
  if ((!is.null(binaryTruth))&(!is.null(inferredGraph))){
    truthVec <- binaryTruth[lower.tri(binaryTruth)]
    inferredGraphVec <- inferredGraph[lower.tri(inferredGraph)]
    FDR <- sum(inferredGraphVec[truthVec==0])/ sum(inferredGraphVec)
    TDR <- sum(inferredGraphVec[truthVec==1])/ sum(truthVec==1)
    xlab2 <- paste('inferred graph, FDR=', round(FDR, digits=3), sep='')
    xlab3 <- paste('true binary graph, TDR=', round(TDR, digits=3), sep='')
    res <- list(FDR=FDR, TDR=TDR)
  }
  if (!is.null(inferredGraph))
    stats::heatmap(inferredGraph, Rowv=NA, Colv=NA, symm=TRUE, col = RColorBrewer::brewer.pal(9,'YlGn'),
            xlab=xlab2)
  if (!is.null(binaryTruth))
    stats::heatmap(binaryTruth, Rowv=NA, Colv=NA, symm=TRUE, col = RColorBrewer::brewer.pal(9,'YlGn'),
            xlab=xlab3)
  return(res)
}

