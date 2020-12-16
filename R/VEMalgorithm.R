#' VEM algorithm to adjust the noisy stochastic block model to an observed dense adjacency matrix
#'
#' \code{fitNSBM()} estimates model parameters of the noisy stochastic block model and provides a clustering of the nodes
#'
#' @details
#' \code{fitNSBM()} supports different probability distributions for the edges and can estimate the number of node blocks
#'
#' @param dataMatrix observed dense adjacency matrix
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
#' @param sbmSize list of parameters determining the size of SBM (the number of latent blocks) to be expored
#' \describe{
#'   \item{\code{Qmin}}{minimum number of latent blocks}
#'   \item{\code{Qmax}}{maximum number of latent blocks}
#'   \item{\code{explor}}{if \code{Qmax} is not provided, then \code{Qmax} is automatically determined as \code{explor} times the number of blocks where the ICL is maximal}
#' }
#' @param filename results are saved in a file with this name (if provided)
#' @param initParam list of parameters that fix the number of initializations
#' \describe{
#'   \item{\code{nbOfTau}}{number of initial points for the node clustering (i. e. for the variational parameters \code{tau})}
#'   \item{\code{nbOfPointsPerTau}}{number of initial points of the latent binary graph}
#'   \item{\code{maxNbOfPasses}}{maximum number of passes through the SBM models, that is, passes from \code{Qmin} to \code{Qmax} or inversely}
#'   \item{\code{minNbOfPasses}}{minimum number of passes through the SBM models}
#' }
#' @param nbCores number of cores used for parallelization
#'
#' @return Returns a list of estimation results for all numbers of latent blocks considered by the algorithm.
#' Every element is a list composed of:
#' \describe{
#'   \item{\code{theta}}{estimated parameters of the noisy stochastic block model; a list with the following elements:
#'   \describe{
#'     \item{\code{pi}}{parameter estimate of pi}
#'     \item{\code{w}}{parameter estimate of w}
#'     \item{\code{nu0}}{parameter estimate of nu0}
#'     \item{\code{nu}}{parameter estimate of nu}
#'    }}
#'   \item{\code{clustering}}{node clustering obtained by the noisy stochastic block model, more precisely, a hard clustering given by the
#'   maximum aposterior estimate of the variational parameters \code{sbmParam$edgeProba}}
#'   \item{\code{sbmParam}}{further results concerning the latent binary stochastic block model. A list with the following elements:
#'   \describe{
#'     \item{\code{Q}}{number of latent blocks in the noisy stochastic block model}
#'     \item{\code{clusterProba}}{soft clustering given by the conditional probabilities of a node to belong to a given latent block.
#'     In other words, these are the variational parameters \code{tau}; (Q x n)-matrix}
#'     \item{\code{edgeProba}}{conditional probabilities \code{rho} of an edges given the node memberships of the interacting nodes; (N_Q x N)-matrix}
#'     \item{\code{ICL}}{value of the ICL criterion at the end of the algorithm}
#'     }}
#'   \item{\code{convergence}}{a list of convergence indicators:
#'   \describe{
#'     \item{\code{J}}{value of the lower bound of the log-likelihood function at the end of the algorithm}
#'     \item{\code{complLogLik}}{value of the complete log-likelihood function at the end of the algorithm}
#'     \item{\code{converged}}{indicates if algorithm has converged}
#'     \item{\code{nbIter}}{number of iterations performed}
#'  }}
#' }
#'
#' @export
#' @examples
#' n <- 10
#' theta <- list(pi= c(0.5,0.5), nu0=c(0,.1),
#'        nu=matrix(c(-2,10,-2, 1,1,1),3,2),  w=c(.5, .9, .3))
#' obs <- rnsbm(n, theta, modelFamily='Gauss')
#' res <- fitNSBM(obs$dataMatrix, sbmSize = list(Qmax=3),
#'        initParam=list(nbOfTau=1, nbOfPointsPerTau=1), nbCores=1)

fitNSBM <- function(dataMatrix, model='Gauss0',
                    sbmSize = list(Qmin=1, Qmax=NULL, explor=1.5),
                    filename=NULL,
                    initParam = list(nbOfTau=NULL, nbOfPointsPerTau=NULL,
                                     maxNbOfPasses=NULL, minNbOfPasses=1),
                    nbCores=parallel::detectCores()){

  directed <- FALSE  # a future version of the R package will include the code for the directed case

  Qmin <- max(c(1, sbmSize$Qmin))
  explor <- if(is.null(sbmSize$explore)) 1.5 else sbmSize$explor
  QmaxIsFlexible <- is.null(sbmSize$Qmax)
  Qmax <- if (QmaxIsFlexible) max(c(4,round(Qmin*explor))) else sbmSize$Qmax

  if (is.null(initParam$nbOfTau)){
    nbOfTau <-  if (model %in%  c('Gauss0EqVar','Gauss0Var1','Gauss01','Gauss2distr','GaussAffil')) 3 else 5
  }else{
    nbOfTau <- initParam$nbOfTau
  }

  if (is.null(initParam$nbOfPointsPerTau)){
    nbOfPointsPerTau <- if (model %in%  c('Gauss0EqVar','Gauss2distr','GaussAffil')) 2 else 6
  }else{
    nbOfPointsPerTau <- initParam$nbOfPointsPerTau
  }
  if (model %in%  c('Gauss0Var1','Gauss01'))
    nbOfPointsPerTau <- 1  # only use Storey to initialize rho

  minNbOfPasses <- max(c(1,initParam$minNbOfPasses))

  maxNbOfPasses <- if (is.null(initParam$maxNbOfPasses)) max(c(10,initParam$minNbOfPasses)) else initParam$maxNbOfPasses

  if (Sys.info()["sysname"]=="Windows"){
    nbCores <- 1
  } else {
    nbCores <- if (nbCores>1) min(c(round(nbCores), parallel::detectCores())) else 1
  }

  doParallelComputing <- (nbCores>1)

  n <- ncol(dataMatrix)
  data <- dataMatrix[lower.tri(diag(n))]
  N <- length(data)
  if (model %in% c('Gauss','Gauss0','Gauss01','GaussEqVar','Gauss0EqVar','Gauss0Var1','Gauss2distr','GaussAffil'))
    modelFamily <- 'Gauss'
  if (model %in% c('Exp','ExpGamma','Gamma'))
    modelFamily <- 'Gamma'

  possibleQvalues <- Qmin:Qmax
  bestSolutionAtQ <- lapply(1:(Qmax-Qmin+1), function(elem) list(convergence=list(J=-Inf))) # list of best solutions for every Q

  currentQindice <- 1
  notConverged <- TRUE
  up <- TRUE  # redundant since up=TRUE iff nbOfPasses is odd
  nbOfPasses <- 1
  while (notConverged){
    Q <- possibleQvalues[currentQindice]
    N_Q <- Q*(Q+1)/2
    cat("-- pass =", nbOfPasses, 'Q =', Q, 'Qmax =', Qmax, '\n')

    ## create list of initial points
    if (nbOfPasses==1){
      ListOfTauRho <- initialPoints(Q, dataMatrix, nbOfTau, nbOfPointsPerTau, modelFamily, model, directed)   # use kmeans, random initializations for gamma parameters, for gamma use more kmeans perturbations??
      if (Q>Qmin){ # add split up solution if Q>Qmin
        additionalTauRho <- initialPointsBySplit(bestSolutionAtQ[[currentQindice-1]]$sbmParam$clusterProba, nbOfTau= round(Q/2), nbOfPointsPerTau, data, modelFamily, model, directed)
        ListOfTauRho$tau <- append(ListOfTauRho$tau, additionalTauRho$tau)
        ListOfTauRho$rho <- append(ListOfTauRho$rho, additionalTauRho$rho)
        if(model=='ExpGamma'){
          ListOfTauRho$nu <- append(ListOfTauRho$nu, additionalTauRho$nu)
        }
      }
    }else{ # nbOfPasses>1
      if (up){ # add split up solution if Q>Qmin
        ListOfTauRho <- initialPointsBySplit(bestSolutionAtQ[[currentQindice-1]]$sbmParam$clusterProba, nbOfTau= max(c(3,round(Q/2))), nbOfPointsPerTau, data, modelFamily, model, directed )
      }else{
        # merge solutions with taudown
        ListOfTauRho <- initialPointsByMerge(bestSolutionAtQ[[currentQindice+1]]$sbmParam$clusterProba, nbOfTau= max(c(5,round(Q/2))), nbOfPointsPerTau, data, modelFamily, model, directed )
      }
    }

    ## launch VEM for those initial points
    M <- length(ListOfTauRho$tau)
    if(doParallelComputing){
      ListOfSolutions <- parallel::mclapply(1:M, function(k){
        mainVEM_Q_par(k, ListOfTauRho, modelFamily, model, data, directed)},
        mc.cores=nbCores)
      ListOfJ <- lapply(ListOfSolutions, function(solutionThisRun) solutionThisRun$convergence$J)
      currentSolution <- ListOfSolutions[[which.max(ListOfJ)]]
      if (currentSolution$convergence$J>bestSolutionAtQ[[currentQindice]]$convergence$J)
        bestSolutionAtQ[[currentQindice]] <- currentSolution
    }else{
      for (s in 1:M){
        cat("-- pass=", nbOfPasses, "Q=",Q,'run=',s,"\n")
        currentInitialPoint <- if (model != 'ExpGamma') list(tau=ListOfTauRho$tau[[s]], rho=ListOfTauRho$rho[[s]]) else
          list(tau=ListOfTauRho$tau[[s]], rho=ListOfTauRho$rho[[s]], nu=ListOfTauRho$nu[[s]])
        currentSolution <- mainVEM_Q(currentInitialPoint, modelFamily, model, data, directed)
        if (currentSolution$convergence$J>bestSolutionAtQ[[currentQindice]]$convergence$J)
          bestSolutionAtQ[[currentQindice]] <- currentSolution
      }
    }

    ## evaluate results of the current pass and decide whether to continue or not
    if (nbOfPasses==1){
      if (Q<Qmax){
        currentQindice <- currentQindice + 1
      }else{ #Q==Qmax
        if (!QmaxIsFlexible){ # fixed Qmax
          if ((Qmax==Qmin)){
            notConverged <- FALSE
          }else{ # fixed Qmax, Qmax>Qmin -> prepare 2nd pass
            bestQ <- getBestQ(bestSolutionAtQ)
            nbOfPasses <- 2
            currentQindice <- currentQindice - 1
            up <- FALSE
          }
        }else{ #nbOfPasses=1, flexible Qmax
          bestQ <- getBestQ(bestSolutionAtQ)  # list $ICL, $Q
          if(bestQ$Q<=round(Qmax/explor)){ # if max is attained for a 'small' Q -> go to 2nd pass
            nbOfPasses <- 2
            currentQindice <- currentQindice - 1
            up <- FALSE
          }else{ # increase Qmax and continue first pass
            oldQmax <- Qmax
            Qmax <- round(bestQ$Q*explor)
            possibleQvalues <- Qmin:Qmax
            bestSolutionAtQ <- append(bestSolutionAtQ, lapply((oldQmax+1):Qmax, function(elem) list(convergence=list(J=-Inf))))
            # bestSolutionAtQ <- append(bestSolutionAtQ, lapply((oldQmax+1):Qmax, function(elem) list(J=-Inf)))
            currentQindice <- currentQindice + 1
          }
        }
      }
    }else{ # nbOfPasses>1
      if((Q>Qmin)&(Q<Qmax)){
        currentQindice <- currentQindice + up - !up
      }else{ # at the end of a pass
        bestQInThisPass <- getBestQ(bestSolutionAtQ)  # list $ICL, $Q
        if(QmaxIsFlexible){
          if(bestQInThisPass$ICL>bestQ$ICL){ # in this pass the solution was improved
            if(bestQInThisPass$Q<=round(Qmax/explor) ){ # the improvement is achieved for a "small" Q -> go to next pass
              if ((Q==Qmin)){ # only stop when Q==Qmin
                bestQ <- bestQInThisPass
                notConverged <- (nbOfPasses<maxNbOfPasses) # nb of max passes attained ?
              }
              nbOfPasses <- nbOfPasses + 1
              currentQindice <- currentQindice - up + !up
              up <- !up
            }else{ # improvement is achieved for a "large" Q -> increase Qmax and continue pass
              oldQmax <- Qmax
              Qmax <- round(bestQ$Q*explor)
              possibleQvalues <- Qmin:Qmax
              bestSolutionAtQ <- append(bestSolutionAtQ, lapply((oldQmax+1):Qmax, function(elem) list(convergence=list(J=-Inf))))
              currentQindice <- currentQindice + 1
              if (!up){ # go to next pass
                bestQ <- bestQInThisPass
                notConverged <- (nbOfPasses<maxNbOfPasses) # nb of max passes attained ?
                nbOfPasses <- nbOfPasses + 1
                up <- TRUE
              }
            }
          }else{ # no improvement -> exit
            notConverged <- (nbOfPasses<minNbOfPasses) #FALSE
          }
        }else{ # fixed Qmax
          if (bestQInThisPass$ICL>bestQ$ICL){ # in this pass the solution was improved -> prepare next pass
            if (Q==Qmin){ # only stop when Q==Qmin
              bestQ <- bestQInThisPass
              notConverged <- (nbOfPasses<maxNbOfPasses) # nb of max passes attained ?
            }
            nbOfPasses <- nbOfPasses + 1
            currentQindice <- currentQindice - up + !up
            up <- !up
          }else{ # no improvement
            if (nbOfPasses>=minNbOfPasses) # min nb of passes attained -> exit
              notConverged <- FALSE
            else{ # continue
              nbOfPasses <- nbOfPasses + 1
              currentQindice <- currentQindice - up + !up
              up <- !up
            }
          }
        }
      }
    }
  }
  if (!is.null(filename)) save(bestSolutionAtQ, file=filename)

  return(bestSolutionAtQ)
}



#' main function of VEM algorithm with fixed number of SBM blocks
#'
#' @param init list of initial points for the algorithm
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
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return list of estimated model parameters and a node clustering; like the output of fitNSBM()
mainVEM_Q <- function(init, modelFamily, model, data, directed){
  nb.iter <- 500
  epsilon <- 1e-6

  Q <- nrow(init$tau)
  N_Q <- if (directed) Q^2 else Q*(Q+1)/2
  n <- ncol(init$tau)

  VE <- list(tau=init$tau, rho=init$rho)
  if (modelFamily=='Gauss'){
    theta.init <- list(w=rep(NA, N_Q),
                       nu0=c(0,1),
                       nu=matrix(c(0,1), N_Q, 2, byrow=TRUE))
  }else{
    if(model=='Exp'){
      theta.init <- list(w=rep(NA, N_Q),
                         nu0=c(1,1),
                         nu=matrix(c(1,1), N_Q, 2, byrow=TRUE))
    }
    if(model=='ExpGamma'){
      theta.init <- list(w=rep(NA, N_Q),
                         nu0=c(1,1),
                         nu=init$nu)
    }
  }
  theta.init$pi <- if (Q>1) rowMeans(VE$tau) else 1
  mstep <- theta.init
  J.old <- -Inf
  Jeval <- list(J=-Inf, complLogLik=-Inf)

  convergence <- list(converged=FALSE, nb.iter.achieved=FALSE)
  it <- 0
  while (sum(convergence==TRUE)==0){
    it <- it+1
    mstep <- Mstep(VE, mstep, model, data, modelFamily, directed)
    VE <- VEstep(VE, mstep, data, modelFamily, directed)
    Jeval <- JEvalMstep(VE, mstep, data, modelFamily, directed)
    convergence$nb.iter.achieved <- (it > nb.iter+1)
    convergence$converged <-  (abs((Jeval$J-J.old)/Jeval$J)< epsilon)
    J.old <- Jeval$J
  }

  solutionThisRun <- list(theta=list(w=mstep$w, nu0=mstep$nu0, nu=mstep$nu, pi=mstep$pi),
                          clustering=NULL,
                          sbmParam=list(Q=Q, clusterProba=VE$tau, edgeProba=VE$rho, ICL=NULL),
                          convergence = list(J=Jeval$J, complLogLik=Jeval$complLogLik,
                                             converged=convergence$converged, nbIter=it))
  solutionThisRun$sbmParam$ICL <- ICL_Q(solutionThisRun, model)
  solutionThisRun$clustering <- if (Q>1) apply(VE$tau,2,which.max) else rep(1, n)

  return(solutionThisRun)
}



#' main function of VEM algorithm for fixed number of latent blocks in parallel computing
#'
#' runs the VEM algorithm the provided initial point
#'
#' @param s indice of initial point in ListOfTauRho to be used for this run
#' @param ListOfTauRho a list of initial points
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
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param directed booelan to indicate whether the model is directed or undirected

#'
#' @return list of estimated model parameters and a node clustering; like the output of fitNSBM()

mainVEM_Q_par <- function(s, ListOfTauRho, modelFamily, model, data, directed){
  currentInitialPoint <- if (model != 'ExpGamma') list(tau=ListOfTauRho$tau[[s]], rho=ListOfTauRho$rho[[s]]) else
    list(tau=ListOfTauRho$tau[[s]], rho=ListOfTauRho$rho[[s]], nu=ListOfTauRho$nu[[s]])
  currentSolution <- mainVEM_Q(currentInitialPoint, modelFamily, model, data, directed)
  return(currentSolution)
}


#' VE-step
#'
#' performs one VE-step, that is, update of tau and rho
#'
#' @param VE list with variational parameters tau and rho
#' @param mstep list with current model parameters and additional auxiliary terms
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param directed booelan to indicate whether the model is directed or undirected
#' @param fix.iter maximal number of iterations for fixed point equation
#'
#' @return updated list \code{VE} with variational parameters tau and rho
VEstep <- function(VE, mstep, data, modelFamily, directed, fix.iter=5){
  Q <- nrow(VE$tau)
  epsilon <- 1e-6
  # compute rho using w
  ind <- 0
  for (q in 1:Q){
    for (l in q:Q){
      ind <- ind + 1
      rho_numerateur <- modelDensity(data, mstep$nu[ind,], modelFamily) * mstep$w[ind]  # taille N
      VE$rho[ind,] <-  rho_numerateur / (rho_numerateur + modelDensity(data, mstep$nu0, modelFamily)*(1-mstep$w[ind]))
    }
  }

  log.w <- log(mstep$w)
  log.w[mstep$w==0] <- 0
  log.1mw <- log(1-mstep$w)
  log.1mw[mstep$w==1] <- 0

  # solve the fixed point equation
  converged <- FALSE
  it <- 0
  while ((!converged)&(it<fix.iter)){
    it <- it+1
    tau.old <- VE$tau
    VE$tau <- tauUpdate(tau.old, log.w, log.1mw, data, VE, mstep, modelFamily, directed)
    converged <- ((max(abs(VE$tau-tau.old))<epsilon))
  }

  return(VE)
}

# piUpdate <- function(VE){
#   Q <- nrow(VE$tau)
#   param_pi <- if (Q>1) rowMeans(VE$tau) else 1
#   return(param_pi)
# }

#' Compute one iteration to solve the fixed point equation in the VE-step
#'
#' @param tau current value of tau
#' @param log.w value of log(w)
#' @param log.1mw value of log(1-w)
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param VE list with variational parameters tau and rho
#' @param mstep list with current model parameters and additional auxiliary terms
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return updated value of tau
tauUpdate <- function (tau, log.w, log.1mw, data, VE, mstep, modelFamily, directed){
  Q <- nrow(tau)
  n <- ncol(tau)

  tau.new <- tau
  if (Q>1){
    logtau <- log(tau)

    for (i in 1:n){
      j.in.ij <- convertNodePair(rep(i,1,n-1), (1:n)[-i], n, directed)	 # n-1
      for (q in 1:Q){
        ind.ql <- convertGroupPair(rep(q,1,Q),1:Q,Q,directed)  # Q
        rho.mat <- VE$rho[ind.ql,j.in.ij]

        sum_j_taujl_rhoijql <- rowSums(tau[,-i]*rho.mat) # Q
        sum_l_taujl_rhoijql <- colSums(tau[,-i]*rho.mat) # n-1

        logf1_l_ijql <- matrix(NA, ncol=n-1, nrow=Q)
        for (indd in 1:Q){
          logf1_l_ijql[indd,] <- log(modelDensity(data[j.in.ij], mstep$nu[ind.ql[indd],], modelFamily))
        }

        logf0 <- log(modelDensity(data[j.in.ij], mstep$nu0, modelFamily))

        term1 <-log.1mw[ind.ql] *rowSums(tau[,-i]) + (log.w[ind.ql]-log.1mw[ind.ql]) *sum_j_taujl_rhoijql
        term1 <- sum(term1)

        term2 <-  logf0 * (colSums(tau[,-i]) - sum_l_taujl_rhoijql)  + colSums(tau[,-i]*rho.mat*logf1_l_ijql)
        term2 <- sum(term2)

        log.rho.mat <- log(rho.mat)
        log.rho.mat[rho.mat==0] <- 0
        log.1mrho.mat <- log(1-rho.mat)
        log.1mrho.mat[rho.mat==1] <- 0
        psi_rho.mat <- rho.mat*log.rho.mat+ (1-rho.mat)*log.1mrho.mat
        term3 <- sum(tau[,-i]*psi_rho.mat)

        logtau[q,i] <- log(mstep$pi[q]) + term1 + term2 - term3
      }

      ## Normalizing in the log space to avoid numerical problems
      ## and going back to exponential with the same normalization
      logtau[,i] <- logtau[,i]-max(logtau[,i])
      tau.new[,i] <- exp(logtau[,i])
      tau.new[,i] <- correctTau(tau.new[,i])
    }
  }
  return(tau.new)
}


#' M-step
#'
#' performs one M-step, that is, update of pi, w, nu, nu0
#'
#' @param VE list with variational parameters tau and rho
#' @param mstep list with current model parameters and additional auxiliary terms
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
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return updated list \code{mstep} with current model parameters and additional auxiliary terms
Mstep <- function(VE, mstep, model, data, modelFamily, directed){
  Q <- nrow(VE$tau)
  n <- ncol(VE$tau)
  N <- length(data)
  N_Q <- nrow(VE$rho)

  mstep$pi <- if (Q>1) rowMeans(VE$tau) else 1

  threshW <- 0.5
  denom_1 <- denom_2 <- denom_3 <- numer_1 <- numer_2 <- numer_3 <- numer_4 <- numer_5 <- numer_6 <- 0
  default_var <- sqrt(log(n))
  ind.ql <- 0
  for (q in 1:Q){
    for (l in q:Q){
      ind.ql <- ind.ql + 1
      tauql <- getTauql(q,l,VE$tau, n, directed)
      tauql_rho <- tauql*VE$rho[ind.ql,]

      s <- sum(tauql_rho)
      mstep$sum_tau[ind.ql] <- sum(tauql)
      mstep$sum_rhotau[ind.ql] <- s

      calcul_w <- s/mstep$sum_tau[ind.ql]
      mstep$w[ind.ql] <- if (calcul_w <=.Machine$double.eps) .Machine$double.eps else calcul_w
      mstep$w[ind.ql] <- if (calcul_w >= 1-.Machine$double.eps) 1-.Machine$double.eps else calcul_w

      if (modelFamily=='Gauss'){
        if (model=='Gauss'){
          if (s > 0){
            mstep$nu[ind.ql,1] <- sum(tauql_rho*data)/s
            sigma <- sqrt(sum(tauql_rho*(data-mstep$nu[ind.ql,1])^2)/s)
            mstep$nu[ind.ql,2] <- if (sigma <.Machine$double.eps) default_var else sigma
          }else{ # s==0
            mstep$nu[ind.ql,] <- c(0,default_var)
          }
          tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
          numer_1 <- numer_1 + sum(tauql_1mrho*data)
          numer_2 <- numer_2 + sum(tauql_1mrho*data^2)
          denom_1 <- denom_1 + sum(tauql_1mrho)
        }

        if(model=='Gauss0'){
          if (s > 0){
            mstep$nu[ind.ql,1] <- sum(tauql_rho*data)/s
            sigma <- sqrt(sum(tauql_rho*(data-mstep$nu[ind.ql,1])^2)/s)
            mstep$nu[ind.ql,2] <- if (sigma <.Machine$double.eps) default_var else sigma
          }else{ # s==0
            mstep$nu[ind.ql,] <- c(0,default_var)
          }
          tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
          numer_2 <- numer_2 + sum(tauql_1mrho*data^2)
          denom_1 <- denom_1 + sum(tauql_1mrho)
        }

        if(model=='Gauss01'){
          if (s > 0){
            mstep$nu[ind.ql,1] <- sum(tauql_rho*data)/s
            sigma <- sqrt(sum(tauql_rho*(data-mstep$nu[ind.ql,1])^2)/s)
            mstep$nu[ind.ql,2] <- if (sigma <.Machine$double.eps) default_var else sigma
          }else{ # s==0
            mstep$nu[ind.ql,] <- c(0,default_var)
          }
          tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
        }

        if (model=='GaussEqVar'){
          mstep$nu[ind.ql,1] <- if (s>0) sum(tauql_rho*data)/s else 0
          numer_4 <- numer_4 + sum(tauql_rho*(data-mstep$nu[ind.ql,1])^2)
          tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
          denom_1 <- denom_1 + sum(tauql_1mrho)
          numer_1 <- numer_1 + sum(tauql_1mrho*data)
          numer_2 <- numer_2 + sum(tauql_1mrho*data^2)
        }

        if (model=='Gauss0EqVar'){
          mstep$nu[ind.ql,1] <- if (s>0) sum(tauql_rho*data)/s else 0
          numer_2 <- numer_2 + sum(tauql_rho*(data-mstep$nu[ind.ql,1])^2)
          tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
          numer_2 <- numer_2 + sum(tauql_1mrho*data^2)
        }

        if (model=='Gauss0Var1'){
          mstep$nu[ind.ql,1] <- if (s>0) sum(tauql_rho*data)/s else 0
        }

        if (model=='Gauss2distr'){
          numer_3 <- numer_3 + sum(tauql_rho*data)  # for mean estimate mu_1
          numer_4 <- numer_4 + sum(tauql_rho*data^2) # for variance estimate sigma_1
          tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
          denom_1 <- denom_1 + sum(tauql_1mrho)
          numer_1 <- numer_1 + sum(tauql_1mrho*data) # for mean estimate mu_0
          numer_2 <- numer_2 + sum(tauql_1mrho*data^2) # for variance estimate sigma_0
        }

        if (model=='GaussAffil'){
          if (q==l){
            numer_3 <- numer_3 + sum(tauql_rho*data)  # for mean estimate mu_in
            numer_4 <- numer_4 + sum(tauql_rho*data^2) # for variance estimate sigma_in
            denom_2 <- denom_2 + sum(tauql_rho)
          }else{
            numer_5 <- numer_5 + sum(tauql_rho*data)  # for mean estimate mu_out
            numer_6 <- numer_6 + sum(tauql_rho*data^2) # for variance estimate sigma_out
            denom_3 <- denom_3 + sum(tauql_rho)
          }
          tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
          denom_1 <- denom_1 + sum(tauql_1mrho)
          numer_1 <- numer_1 + sum(tauql_1mrho*data) # for mean estimate mu_0
          numer_2 <- numer_2 + sum(tauql_1mrho*data^2) # for variance estimate sigma_0
        }
      }

      if (modelFamily=='Gamma'){
        if (s > 0){
          if (model=='Exp'){
            mstep$nu[ind.ql,][2] <- s/sum(tauql_rho*data)
          }
          if (model=='ExpGamma'){
            L_ql <- sum(tauql_rho*log(data))/s
            M_ql <- sum(tauql_rho*data)/s
            mstep$nu[ind.ql,] <- emv_gamma(L_ql, M_ql, mstep$nu[ind.ql,])
          }
        }else{ # if no observations: use random values
          mstep$nu[ind.ql,][2] <- stats::rgamma(1,2.5,0.5)
        }
        tauql_1mrho <- tauql*(1-VE$rho[ind.ql,])
        numer_2 <- numer_2 + sum(tauql_1mrho)
        denom_1 <- denom_1 + sum(tauql_1mrho*data)
      }

    } # end for l
  }# end for q
  if (modelFamily=='Gauss'){
    if (model=='GaussEqVar'){
      if (denom_1>0){
        mstep$nu0[1] <- numer_1/denom_1
      }
      sigma <- if (denom_1>0) sqrt((denom_1*(numer_2/denom_1-mstep$nu0[1]^2) + numer_4)/N) else default_var
      mstep$nu0[2] <- sigma
      mstep$nu[,2] <- sigma
    }

    if (model=='Gauss0EqVar'){
      mstep$nu0[2] <- sqrt(numer_2/N)
      mstep$nu[,2] <- mstep$nu0[2]
    }

    if (model=='Gauss2distr'){
      if (denom_1>0){
        mstep$nu0[1] <- numer_1/denom_1
        mstep$nu0[2] <- sqrt(numer_2/denom_1 - mstep$nu0[1]^2)
      }else{
        mstep$nu0 <- c(0,default_var)
      }
      if (denom_1<N){
        mstep$nu[,1] <- numer_3/(N-denom_1)
        mstep$nu[,2] <- sqrt(numer_4/(N-denom_1) - mstep$nu[1,1]^2)
      }else{
        mstep$nu <- matrix(c(0,default_var),N_Q,2,byrow=TRUE)
      }
    }
    if (model=='Gauss0'){
      mstep$nu0[2] <- if (denom_1>0) sqrt(numer_2/denom_1) else default_var
    }

    if (model=='Gauss'){
      if (denom_1>0){
        mstep$nu0[1] <- numer_1/denom_1
        mstep$nu0[2] <- sqrt(numer_2/denom_1 - mstep$nu0[1]^2)
      }else{
        mstep$nu0 <- c(0,default_var)
      }
    }

    if (model=='GaussAffil'){
      if (denom_1>0){ # nu_0
        mstep$nu0[1] <- numer_1/denom_1
        mstep$nu0[2] <- sqrt(numer_2/denom_1 - mstep$nu0[1]^2)
      }else{
        mstep$nu0 <- c(0,default_var)
      }
      ind_qq_in <- convertGroupPair(1:Q,1:Q,Q,directed)
      if (denom_2>0){ # nu_in
        mstep$nu[ind_qq_in,1] <- numer_3/denom_2
        mstep$nu[ind_qq_in,2] <- sqrt(numer_4/denom_2 - mstep$nu[1,1]^2)
      }else{
        mstep$nu[ind_qq_in,] <- matrix(c(0,default_var),Q,2,byrow=TRUE)
      }
      if (Q>1){
        ind_qq_out <- setdiff(1:N_Q, ind_qq_in)
        if (denom_3>0){ # nu_out
          mstep$nu[ind_qq_out,1] <- numer_5/denom_3
          mstep$nu[ind_qq_out,2] <- sqrt(numer_6/denom_3 - mstep$nu[2,1]^2)
        }else{
          mstep$nu[ind_qq_out,] <- matrix(c(0,default_var),Q,2,byrow=TRUE)
        }
      }
    }

    if ((mstep$nu0[2]<1e-7)|(is.na(mstep$nu0[2])) )
    {print(paste("nu0=", mstep$nu0[2])) ;
      mstep$nu0[2] <- default_var
    }
    varZero <- ((mstep$nu[,2] < 1e-7) | (is.na(mstep$nu[,2])))
    if (sum(varZero) > 0)
      mstep$nu[varZero,2] <- default_var
  }

  if (modelFamily=='Gamma'){
    mstep$nu0[2] <- if (denom_1>0) numer_2/denom_1 else 1/mean(data[data<stats::median(data)])
  }
  return(mstep)
}



#' Perform one iteration of the Newton-Raphson to compute the MLE of the parameters of the Gamma distribution
#'
#' @param param current parameters of the Gamma distribution
#' @param L weighted mean of log(data)
#' @param M weighted mean of the data
#'
#' @return updated parameters of the Gamma distribution
update_newton_gamma <- function(param, L, M){
  gradient <- c(log(param[2])-digamma(param[1])+L, param[1]/param[2]-M)
  w <- trigamma(param[1])
  detH <- 1/(1-param[1]*w)
  a.new <- min(c(param[1] - detH*sum(param*gradient), 150))
  b.new <- param[2] - detH*sum(c(param[2],param[2]^2*w)*gradient)
  return(c(a.new,b.new))
}

#' evaluate the objective in the Gamma model
#'
#' @param param parameters of the Gamma distribution
#' @param L weighted mean of log(data)
#' @param M weighted mean of the data
#'
#' @return value of the lower bound of the log-likelihood function
J.gamma <- function(param, L, M){
  return(param[1]*log(param[2]) - log(gamma(param[1])) + (param[1]-1)*L - param[2]*M)
}


#' compute the MLE in the Gamma model using the Newton-Raphson method
#'
#' @param L weighted mean of log(data)
#' @param M weighted mean of the data
#' @param param.old parameters of the Gamma distribution
#' @param epsilon threshold for the stopping criterion
#' @param nb.iter.max maximum number of iterations
#'
#' @return updated parameters of the Gamma distribution
emv_gamma <- function(L, M, param.old, epsilon=1e-3, nb.iter.max=10){
  logvrais.old <- J.gamma(param.old, L, M)
  notConverged <- TRUE
  nb.iter <- 0
  while ((notConverged) & (nb.iter <= nb.iter.max)){
    param.new <- update_newton_gamma(param.old, L, M)
    if (sum(param.new>0)==2){
      # check convergence criterion
      logvrais.new <- J.gamma(param.new, L, M)
      notConverged <- (abs((logvrais.new-logvrais.old)/logvrais.new)>epsilon)
      param.old <- param.new
      logvrais.old <- logvrais.new
    }
    else{
      notConverged <- FALSE
    }
    nb.iter <- nb.iter + 1
  }
  return(param.old)
}


#' evaluation of the objective in the Gauss model
#'
#' @param VE list with variational parameters tau and rho
#' @param mstep list with current model parameters and additional auxiliary terms
#' @param data data vector in the undirected model, data matrix in the directed model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param directed booelan to indicate whether the model is directed or undirected
#'
#' @return value of the ELBO and the complete log likelihood function
JEvalMstep  <- function(VE, mstep, data, modelFamily, directed){
  N_Q <- nrow(VE$rho)
  Q <- nrow(VE$tau)
  n <- ncol(VE$tau)

  log.tau <- log(VE$tau)
  log.tau[VE$tau==0] <- 0
  log.pi  <- log(mstep$pi)
  log.pi[mstep$pi==0]   <- 0

  log.rho <- log(VE$rho)
  log.rho[VE$rho==0] <- 0
  log.1mrho <- log(1-VE$rho)
  log.1mrho[VE$rho==1] <- 0

  log.w <- log(mstep$w)
  log.w[mstep$w==0] <- 0
  log.1mw <- log(1-mstep$w)
  log.1mw[mstep$w==1] <- 0

  term1 <- sum(log.1mw*mstep$sum_tau + (log.w-log.1mw)*mstep$sum_rhotau)
  term2 <- term3 <- rep(NA, N_Q)

  ind.ql <- 0
  for (q in 1:Q){
    for (l in q:Q){
      ind.ql <- ind.ql+1
      tauql <- getTauql(q,l, VE$tau, n, directed)

      logf1 <- log(modelDensity(data, mstep$nu[ind.ql,], modelFamily))
      logf0 <- log(modelDensity(data, mstep$nu0, modelFamily))

      term2[ind.ql] <- sum(tauql*VE$rho[ind.ql,]*logf1 + tauql*(1-VE$rho[ind.ql,])*logf0)

      term3[ind.ql] <- sum(tauql*(VE$rho[ind.ql,] * log.rho[ind.ql,] + (1-VE$rho[ind.ql,])*log.1mrho[ind.ql,]))
    }
  }
  term2_sum <- sum(term2)
  term3_sum <- sum(term3)

  J <- sum( log.pi%*% VE$tau )- sum( VE$tau * log.tau ) + term1 + term2_sum - term3_sum
  complLogLik <- sum( log.pi%*% VE$tau )+ term1 + term2_sum

  return(list(J=J, complLogLik=complLogLik))
}



#' computation of the Integrated Classification Likelihood criterion
#'
#' computation of the Integrated Classification Likelihood criterion for a result provided by mainVEM_Q()
#'
#' @param solutionThisRun result provided by mainVEM_Q()
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
#'
#' @return value of the ICL criterion
ICL_Q <- function(solutionThisRun, model){
  N_Q <- nrow(solutionThisRun$sbmParam$edgeProba)
  Q <- solutionThisRun$sbmParam$Q
  n <- ncol(solutionThisRun$sbmParam$clusterProba)
  N <- ncol(solutionThisRun$sbmParam$edgeProba)
  dimW <- N_Q
  if (model %in% c('Gauss0', 'ExpGamma')){
    dimH0 <- 1
    dimH1 <- 2*N_Q
  }
  if (model=='Gauss'){
    dimH0 <- 2
    dimH1 <- 2*N_Q
  }
  if (model=='Gauss01'){
    dimH0 <- 0
    dimH1 <- 2*N_Q
  }
  if (model=='Gauss0Var1'){
    dimH0 <- 0
    dimH1 <- N_Q
  }
  if (model=='GaussEqVar'){
    dimH0 <- 2 # = 1 gaussian mean under H0 + 1 common variance for all gaussians in the model
    dimH1 <- N_Q # nb of gaussian means under H1
  }
  if (model=='Gauss0EqVar'){
    dimH0 <- 1 # 1 common variance for all gaussians in the model
    dimH1 <- N_Q # nb of gaussian means under H1
  }
  if (model=='GaussAffil'){
    dimH0 <- 2
    dimH1 <- 4
  }
  if (model=='Gauss2distr'){
    dimH0 <- 2
    dimH1 <- 2
  }
  if (model=='Exp'){
    dimH0 <- 1
    dimH1 <- N_Q
  }
  dimParam <- dimW + dimH0 + dimH1

  penalty <- (Q-1)*log(n)+ dimParam*log(N)
  ICL <- 2*solutionThisRun$convergence$complLogLik - penalty
  return(ICL)
}



#' optimal number of SBM blocks
#'
#' returns the number of SBM blocks that maximizes the ICL
#'
#' @param bestSolutionAtQ output of \code{fitNSBM()}, i.e. a list of estimation results for varying number of latent blocks
#'
#' @return a list the maximal value of the ICL criterion among the provided solutions along with the best number of latent blocks
#' @export
#'
#' @examples
#' # res_gauss is the output of a call of fitNSBM()
#' getBestQ(res_gauss)

getBestQ <- function(bestSolutionAtQ){
  L <- length(bestSolutionAtQ)
  Qmin <- bestSolutionAtQ[[1]]$sbmParam$Q
  Qmax <- bestSolutionAtQ[[L]]$sbmParam$Q
  possibleQvalues <- Qmin:Qmax
  listICL <- sapply(bestSolutionAtQ, function(elem) elem$sbmParam$ICL)
  return(list(Q=possibleQvalues[which.max(listICL)], ICL = max(listICL)))
}

#' plot ICL curve
#'
#' @param res output of fitNSBM()
#'
#' @return figure of ICL curve
#' @export
#' @examples
#' # res_gauss is the output of a call of fitNSBM()
#' plotICL(res_gauss)
plotICL <- function(res){
  Q <- ICL <- NULL
  hatQ <- getBestQ(res)
  dfICL <- data.frame(ICL = sapply(res, function(el) el$sbmParam$ICL),
                      Q = sapply(res, function(el) el$sbmParam$Q))
  g <- ggplot2::ggplot(dfICL, ggplot2::aes(x=Q, y=ICL)) +
    ggplot2::geom_line() +
    ggplot2::xlab(paste('best number of blocks:', hatQ))
  return(g)
}

