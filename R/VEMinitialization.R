
#' compute a list of initial points for the VEM algorithm
#'
#' compute a list of initial points of tau and rhofor the VEM algorithm
#' for a given number of blocks; returns nbOfTau*nbOfPointsPerTau inital points
#'
#' @param Q number of latent blocks in the noisy stochastic block model
#' @param dataMatrix observed dense adjacency matrix
#' @param nbOfTau number of initializations for the latent block memberships
#' @param nbOfPointsPerTau number of initializations of the latent binary graph
#' associated with each initial latent block memberships
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return list of inital points of tau and rho of length nbOfTau*nbOfPointsPerTau
initialPoints <- function(Q, dataMatrix, nbOfTau, nbOfPointsPerTau, modelFamily, model,
                          directed){
  n <- nrow(dataMatrix)
  tau <- initialTau(Q, dataMatrix, nbOfTau, percentageOfPerturbation=0.25, directed)
  data <- if (directed) dataMatrix else dataMatrix[lower.tri(diag(n))]
  init <- initialRho(tau, nbOfPointsPerTau, data, modelFamily, model, directed)
  return(init)
}


#' compute intial values for tau
#'
#' returns a list of length nbOfTau of initial points for tau using spectral clustering with absolute values, kmeans and random perturbations of these points
#'
#' @param Q number of latent blocks in the noisy stochastic block model
#' @param dataMatrix observed dense adjacency matrix
#' @param nbOfTau number of initializations for the latent block memberships
#' @param percentageOfPerturbation percentage of node labels that are perturbed to obtain further initial points
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return a list of length nbOfTau of initial points for tau
initialTau <- function(Q, dataMatrix, nbOfTau, percentageOfPerturbation, directed){
  n <- nrow(dataMatrix)

  if (Q==1){
    listOfTau <- lapply(1:nbOfTau, function(x) (matrix(1,1,n)) )
  }else{
    listOfTau <- vector("list", nbOfTau)

    k <- 1
    # absolute spectral clustering
    if (!directed){
      trySpect <- try(spectralClustering(dataMatrix, Q), silent=TRUE)
    }else{ # symetrize the adjacency matrix
      A <- dataMatrix + t(dataMatrix)
      trySpect <- try(spectralClustering(A, Q), silent=TRUE)
    }
    if (is.numeric(trySpect)){
      listOfTau[[k]] <- apply(t(classInd(trySpect,Q)),2,correctTau)
      k <- k + 1
    }

    # k-means
    tryKmeans <- try(stats::kmeans(dataMatrix,Q,nstart=50),silent=TRUE)
    if (is.numeric(tryKmeans$cluster)){
      listOfTau[[k]] <-  apply(t(classInd(tryKmeans$cluster,Q)),2,correctTau)
      k <- k + 1
    }

    # if directed : also apply kmeans on the rows of the adjacency matrix
    if (directed){
      tryKmeans <- try(stats::kmeans(t(dataMatrix),Q,nstart=50),silent=TRUE)
      if (is.numeric(tryKmeans$cluster)){
        listOfTau[[k]] <- apply(t(classInd(tryKmeans$cluster,Q)),2,correctTau)
        k <- k + 1
      }
    }

    # perturbation of spectral clustering and the first k-means solution by replacing some labels at random
    nbswitch <- min(c(max(c(round(n*percentageOfPerturbation),2)),n))
    s <- k
    while(k <= nbOfTau){
      listOfTau[[k]] <- listOfTau[[k%%s + 1]]
      label2switch <- sample(1:n, nbswitch)
      newClasses <- sample(1:Q, nbswitch, replace=TRUE)
      listOfTau[[k]][ , label2switch] <- t(classInd(newClasses, Q))
      listOfTau[[k]] <- apply(listOfTau[[k]], 2, correctTau)
      k <- k + 1
    }
  }
  return(listOfTau)
}




#' spectral clustering with absolute values
#'
#' performs absolute spectral clustering of an adjacency matrix
#'
#' @param A adjacency matrix
#' @param K number of desired clusters
#'
#' @return a vector containing a node clustering into K groups
spectralClustering <- function(A, K){
  n <- nrow(A)
  matrixDmin1_2 <- diag(1/sqrt(rowSums(abs(A)))) # use absolute value of the adjacency matrix to avoid numrical issues
  Labs <- matrixDmin1_2 %*% A %*%matrixDmin1_2
  specabs <- eigen(Labs)
  index <- order(abs(specabs$values),decreasing = FALSE)[(n-K+1):n]
  U <- specabs$vectors[,index]
  clustering <- stats::kmeans(U,K,nstart=100)$cluster
  return(clustering)
}

#' compute initial values of rho
#'
#' for every provided initial point of tau nbOfPointsPerTau initial values of rho are computed
#' in the Gamma model also initial values of nu are computed
#'
#' @param listOfTau output of initialTau()
#' @param nbOfPointsPerTau number of initializations of the latent binary graph associated with each initial latent block memberships
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return list of inital points of tau and rho
initialRho <- function(listOfTau, nbOfPointsPerTau, data, modelFamily, model, directed){
  N <- length(data)
  Q <- nrow(listOfTau[[1]])
  n <- ncol(listOfTau[[1]])
  N_Q <- if (directed) Q^2 else Q*(Q+1)/2

  nbOfTau <- length(listOfTau)
  M <- nbOfTau*nbOfPointsPerTau
  init <- list(tau=lapply(1:M, function(k) listOfTau[[ceiling(k/nbOfPointsPerTau)]]),
               rho=lapply(1:M, function(k) matrix(1, N_Q, N))  # for every tau keep one initialization with all rho equal to 1
  )

  if (modelFamily=='Gauss'){
    nu0 <- c(0,1)
    pvalues <- 2*(1-stats::pnorm(abs(data)))
    for (k in 1:M){
      if (k%%nbOfPointsPerTau == 0){ # initialisation with Storey
        nu <- matrix(c(0,1), N_Q, 2, byrow=TRUE)
        w <- rep(0, N_Q)
        ind_ql <- 0
        for (q in 1:Q){
          for (l in q:Q){
            ind_ql <- ind_ql + 1
            tauql <- getTauql(q, l, init$tau[[k]], n, directed)

            ss <- sum(round(tauql))
            if (ss>0)
              w[ind_ql] <- max(c(1-2*sum(round(tauql)*(pvalues>0.5))/ss,0))
            tauql <- tauql*(pvalues<=0.5)
            s <- sum(tauql)
            if (s>0)
              nu[ind_ql,1] <- sum(tauql*data)/s
          }
        }
        init$rho[[k]] <- getRho(Q, w, nu0, nu, data, modelFamily)
      }else{
        if (k%%nbOfPointsPerTau > 1){ # initialisation at random
          w <- stats::rbeta(N_Q,2,2)
          init$rho[[k]] <- t(sapply(1:N_Q, function(ind_ql) stats::rbeta(N, 1, (1-w[ind_ql])/w[ind_ql])))
        }
      }
    }
  }

  if (modelFamily=='Gamma'){
    dataMean <- mean(data)
    if (model=='Exp'){
      for (k in 1:M){
        if (k%%nbOfPointsPerTau != 0){
          w <-  stats::rbeta(N_Q,2,2)
          nu0 <- c(1, stats::rexp(1, dataMean))
          nu <- matrix( c(rep(1,N_Q), stats::rgamma(N_Q, 1/dataMean, 1)), N_Q, 2)
          init$rho[[k]] <- getRho(Q, w, nu0, nu, data, modelFamily)
        }
      }
    }else{ # 'ExpGamma'
      init$nu <- lapply(1:M, function(k) matrix(1, N_Q, 2))
      dataVar <- if (directed) (n^2-1)/N*stats::var(data) else stats::var(lower.tri(data))
      for (k in 1:M){
        if (k%%nbOfPointsPerTau != 0){ # else keep all rho=1
          w <- stats::rbeta(N_Q,2,2)
          nu0 <- c(1, stats::rexp(1, dataMean))
          m <- rep(-1, N_Q)
          while(sum(m<0)>0)
            m[m<0] <- stats::rnorm(sum(m<0), dataMean, sqrt(dataMean))
          v <- rep(-1, N_Q)
          while(sum(v<0)>0)
            v[v<0] <- stats::rnorm(sum(v<0), dataVar, sqrt(dataMean))
          init$nu[[k]][,1] <- m^2/v
          init$nu[[k]][,2] <- m/v
          init$rho[[k]] <- getRho(Q, w, nu0, init$nu[[k]], data, modelFamily)
        }
      }
    }
  }
  return(init)
}



#' compute rho associated with given values of w, nu0 and nu
#'
#' @param Q number of latent blocks in the noisy stochastic block model
#' @param w weight parameter in the noisy stochastic block model
#' @param nu0 null parameter in the noisy stochastic block model
#' @param nu alternative parameter in the noisy stochastic block model
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#'
#' @return a matrix of conditional probabilities of an edge given the node memberships of the interacting nodes
getRho <- function(Q, w, nu0, nu, data, modelFamily){
  N_Q <- length(w)
  N <- length(data)
  rho <-  matrix(NA, nrow=N_Q, ncol=N)
  ind_ql <- 0
  for (q in 1:Q){
    for (l in q:Q){
      ind_ql <- ind_ql + 1
      rhoNumerator <- modelDensity(data, nu[ind_ql,], modelFamily) * w[ind_ql]  # taille N
      rho[ind_ql,] <- rhoNumerator / (rhoNumerator + modelDensity(data, nu0, modelFamily)*(1-w[ind_ql]))   # taille N
    }
  }
  return(rho)
}




#' convert a clustering into a 0-1-matrix
#'
#' @param cl cluster in vector form
#' @param nbClusters number of clusters
#'
#' @return a 0-1-matrix encoding the clustering
classInd <- function (cl, nbClusters){
  nbClusters <- max(nbClusters, max(cl))
  n <- length(cl)
  x <- matrix(0, n, nbClusters)
  x[(1:n) + n*(cl-1)] <- 1
  return(x)
}


#' Construct initial values with Q groups by splitting groups of a solution obtained with Q-1 groups
#'
#' @param tau_Qm1 tau for a model with Q-1 latent blocks
#' @param nbOfTau number of initializations for the latent block memberships
#' @param nbOfPointsPerTau number of initializations of the latent binary graph associated with each initial latent block memberships
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return list of inital points of tau and rho of length nbOfTau*nbOfPointsPerTau
initialPointsBySplit <- function(tau_Qm1, nbOfTau, nbOfPointsPerTau, data, modelFamily, model, directed){
  tau <- tauUp(tau_Qm1, nbOfTau)
  init <- initialRho(tau, nbOfPointsPerTau, data, modelFamily, model, directed)
  return(init)
}

#' Create new values of tau by splitting groups of provided tau
#'
#' Create nbOfSplits (or all) new values of tau by splitting nbOfSplits (or all) groups of provided tau
#'
#' @param tau soft node clustering
#' @param nbOfSplits number of required splits of blocks
#'
#' @return a list of length nbOfSplits (at most) of initial points for tau
tauUp <- function(tau, nbOfSplits=1){
  n <- ncol(tau)
  Qold <- nrow(tau) # value of Q at previous step (for which a solution is available)

  if (Qold==1){
    listOfTau <- lapply(1:nbOfSplits, function(k) t(gtools::rdirichlet(n, c(0.7,0.7))) )
  }else{ # Qold>=2
    if (nbOfSplits>=Qold){ # split all the groups once
      listOfTau <- lapply(1:Qold, function(q) addRowToTau(tau,q))
    }else{ # always split the largest group
      largestGroup <- which.max(apply(tau,1,sum))
      if (nbOfSplits>=2){ # and split the component with largest entropy
        entropy <- tau*log(tau)
        entropy[is.na(entropy)] <- 0
        smallestEntropy <- which.min(apply(entropy,1,sum))
      }
      groupsToSplit <- if (nbOfSplits>2) sample((1:Qold)[-c(largestGroup,smallestEntropy)],nbOfSplits-2,replace=F) else NULL

      listOfTau <- lapply(c(largestGroup, smallestEntropy, groupsToSplit), function(q) addRowToTau(tau,q))
    }
  }

  return(listOfTau)
}




#' split group q of provided tau randomly into two into
#'
#' @param tau provided tau
#' @param q indice of group to split
#'
#' @return new tau
addRowToTau <- function(tau, q){
  n <- ncol(tau)

  newTau <- tau
  newTau[q,] <- newTau[q,]*stats::runif(n)
  newTau <- rbind(newTau,1-newTau[q,])
  newTau <- apply(newTau, 2, function(col) col/sum(col))
  return(newTau)
}

#' Construct initial values with Q groups by meging groups of a solution obtained with Q+1 groups
#'
#' @param tau_Qp1 tau for a model with Q+1 latent blocks
#' @param nbOfTau number of initializations for the latent block memberships
#' @param nbOfPointsPerTau number of initializations of the latent binary graph associated with each initial latent block memberships
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return list of inital points of tau and rho of length nbOfTau*nbOfPointsPerTau
initialPointsByMerge <- function(tau_Qp1, nbOfTau, nbOfPointsPerTau, data, modelFamily, model, directed){
  tau <- tauDown(tau_Qp1, nbOfTau)
  init <- initialRho(tau, nbOfPointsPerTau, data, modelFamily, model, directed)
  return(init)
}

#' Create new initial values by merging pairs of groups of provided tau
#'
#' Create nbOfMerges new initial values by merging nbOfMerges (or all possible) pairs of groups of provided tau
#'
#' @param tau soft node clustering
#' @param nbOfMerges number of required merges of blocks
#'
#' @return a list of length nbOfMerges (at most) of initial points for tau
tauDown <- function(tau, nbOfMerges){
  n <- ncol(tau)
  Qold <- nrow(tau)
  if (Qold==2){
    listOfTau <- list(matrix(1,1,n))
  } else {
    listOfTau <- list()
    if (nbOfMerges>=Qold*(Qold-1)/2){ # merge all the possible pairs of groups (q,l) with q<l
      for (q in 1:(Qold-1)){
        for (l in ((q+1):Qold)){
          newTau <- tau
          newTau <- rbind(newTau[-c(q,l),], newTau[q,]+newTau[l,])
          listOfTau <- append(listOfTau, list(newTau))
        }
      }
    } else {
      # start by merging the 2 components with largest entropies
      entropy <- tau*log(tau)
      entropy[is.na(entropy)] <- 0
      entropyPerGroup <- apply(entropy,1,sum)
      q <- order(entropyPerGroup)[1]
      l <- order(entropyPerGroup)[2]
      newTau <- tau
      newTau <- rbind(newTau[-c(q,l),], newTau[q,]+newTau[l,])
      listOfTau <- append(listOfTau, list(newTau))
      # merge a number nbOfMerges-1 of pairs of groups (q,l) q<l which are randomly chosen
      groupsToMerge <- sample(1:(Qold*(Qold-1)/2), nbOfMerges-1, replace=FALSE)
      for (ind_ql in groupsToMerge){
        newTau <- tau
        q <- convertGroupPairIdentifier(ind_ql, Qold)[1]
        l <- convertGroupPairIdentifier(ind_ql, Qold)[2]
        newTau <- rbind(newTau[-c(q,l),], newTau[q,]+newTau[l,])
        listOfTau <- append(listOfTau, list(newTau))
      }
    } # end if/else nbOfMerges>=Q*(Q-1)/2
  } # end if/else Q==2
  return(listOfTau)
}




