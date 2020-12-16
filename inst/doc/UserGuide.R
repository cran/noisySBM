## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)

## -----------------------------------------------------------------------------
library(noisySBM)

## -----------------------------------------------------------------------------
directed <- FALSE

## -----------------------------------------------------------------------------
theta <- list(pi=NULL, w=NULL, nu0=NULL, nu=NULL)

## -----------------------------------------------------------------------------
 Q <- 2

## -----------------------------------------------------------------------------
theta$pi <- c(2,1)/3

## -----------------------------------------------------------------------------
theta$w <- c(.8,.1,.9)

## -----------------------------------------------------------------------------
theta$nu0 <- c(0,1)

## -----------------------------------------------------------------------------
theta$nu <- matrix(c(-1,5,-1, 1,1,1), ncol=2)
theta$nu

## -----------------------------------------------------------------------------
obs <- rnsbm(n = 10, theta, modelFamily = 'Gauss')
round(obs$dataMatrix, digits = 2)

## -----------------------------------------------------------------------------
obs$latentAdj

## -----------------------------------------------------------------------------
obs$latentZ

## -----------------------------------------------------------------------------
plotGraphs(obs$dataMatrix, binaryTruth=obs$latentAdj)

## -----------------------------------------------------------------------------
theta$w <- c(.8,.1,.2,.9)
theta$nu <- matrix(c(-1,5,10,-1, 1,1,1,1), ncol=2)
theta$nu

## -----------------------------------------------------------------------------
obs <- rnsbm(n = 10, theta, modelFamily = 'Gauss', directed = TRUE)
plotGraphs(obs$dataMatrix, binaryTruth=obs$latentAdj)
round(obs$dataMatrix, digits = 2)
obs$latentAdj

## -----------------------------------------------------------------------------
theta$pi <- rep(1/3, 3)

## -----------------------------------------------------------------------------
theta$w <- c(.1, .8, .8, .1, .8, .1)

## -----------------------------------------------------------------------------
theta$nu0 <- c(1,2)

## -----------------------------------------------------------------------------
theta$nu <- matrix(c(1, 1, 1, 1, 1, 1,  1, .5, .5, 1, .5, 1), ncol=2)
theta$nu

## -----------------------------------------------------------------------------
obs <- rnsbm(n = 10, theta, modelFamily = 'Gamma')
plotGraphs(obs$dataMatrix, binaryTruth=obs$latentAdj)
round(obs$dataMatrix, digits = 2)
obs$latentAdj
obs$latentZ

## -----------------------------------------------------------------------------
theta <- list(pi=1, w=.5, nu0=1, nu=5)

## -----------------------------------------------------------------------------
obs <- rnsbm(n = 10, theta, modelFamily = 'Poisson', directed = TRUE)
plotGraphs(obs$dataMatrix, binaryTruth=obs$latentAdj)
round(obs$dataMatrix)
obs$latentAdj
obs$latentZ

## -----------------------------------------------------------------------------
set.seed(1)
theta1 <- list(pi=c(.5,.5), w=c(.8,.1,.2), nu0=c(0,1), nu=matrix(c(-1,5,10, 1,1,1), ncol=2))
obs1 <- rnsbm(n=30, theta1)

## ---- eval=FALSE--------------------------------------------------------------
#  res_gauss <- fitNSBM(obs1$dataMatrix, nbCores=2)

## -----------------------------------------------------------------------------
res_gauss[[2]]$theta

## -----------------------------------------------------------------------------
res_gauss[[2]]$clustering

## -----------------------------------------------------------------------------
res_gauss[[2]]$sbmParam$Q

## -----------------------------------------------------------------------------
res_gauss[[2]]$sbmParam$clusterProba

## -----------------------------------------------------------------------------
res_gauss[[2]]$sbmParam$edgeProba[,1:5]

## -----------------------------------------------------------------------------
dim(res_gauss[[2]]$sbmParam$edgeProba)

## -----------------------------------------------------------------------------
res_gauss[[2]]$sbmParam$ICL

## -----------------------------------------------------------------------------
res_gauss[[2]]$convergence

## -----------------------------------------------------------------------------
theta2 <- list(pi=rep(1/2, 2), 
               w=c(.2, .8, .2), 
               nu0=c(1,5), 
               nu=matrix(c(1, 1, 1,   1, 1/3, 1), ncol=2))
set.seed(2)
obs2 <- rnsbm(n=60, theta2, modelFamily='Gamma')

## ---- eval=FALSE--------------------------------------------------------------
#  res_exp <- fitNSBM(obs2$dataMatrix, model='Exp', nbCores = 2)
#  
#  set.seed(3)
#  res_gamma <- fitNSBM(obs2$dataMatrix, model='ExpGamma', nbCores = 2)

## -----------------------------------------------------------------------------
Q <- 2
res_exp[[Q]]$theta
res_gamma[[Q]]$theta

## -----------------------------------------------------------------------------
getBestQ(res_gauss)

## -----------------------------------------------------------------------------
plotICL(res_gauss)

## -----------------------------------------------------------------------------
getBestQ(res_exp)
getBestQ(res_gamma)

## ---- warning=FALSE-----------------------------------------------------------
dfICLexp <- data.frame(ICL = sapply(res_exp, function(el) el$sbmParam$ICL), 
                    Q = sapply(res_exp, function(el) el$sbmParam$Q), model='Exp')
dfICLgam <- data.frame(ICL = sapply(res_gamma, function(el) el$sbmParam$ICL), 
                    Q = sapply(res_gamma, function(el) el$sbmParam$Q), model='ExpGamma')
dfICL <- rbind(dfICLexp, dfICLgam)
ggplot(dfICL, aes(x=Q, y=ICL, group=model, colour=model) ) +
  geom_line() +
  xlim(1,4)

## ---- eval=FALSE--------------------------------------------------------------
#  res4 <- fitNSBM(obs1$dataMatrix, sbmSize=list(Qmin=2, Qmax=2), nbCores=2)

## -----------------------------------------------------------------------------
ARI(res_gauss[[3]]$clustering, obs1$latentZ)

## -----------------------------------------------------------------------------
ARI(res_gauss[[2]]$clustering, res_gauss[[3]]$clustering)

## -----------------------------------------------------------------------------
ARI(res_exp[[2]]$clustering, obs2$latentZ)
ARI(res_gamma[[2]]$clustering, obs2$latentZ)

## -----------------------------------------------------------------------------
resGraph <- graphInference(obs1$dataMatrix, res_gauss[[2]]$clustering, res_gauss[[2]]$theta, alpha = .05, modelFamily='Gauss')
resGraph$A[1:5, 1:5]

## -----------------------------------------------------------------------------
resGraph$qval[1:10]

## -----------------------------------------------------------------------------
length(resGraph$qval)

## -----------------------------------------------------------------------------
plotGraphs(obs1$dataMatrix, resGraph$A, obs1$latentAdj)

## -----------------------------------------------------------------------------
resGraph2 <- graphInference(obs2$dataMatrix, res_exp[[2]]$clustering, res_exp[[2]]$theta, alpha = .05, modelFamily='Gamma')
resGraph2$A[1:5, 1:5]

## -----------------------------------------------------------------------------
resGraph2$qval[1:10]

## -----------------------------------------------------------------------------
length(resGraph2$qval)

## -----------------------------------------------------------------------------
plotGraphs(obs2$dataMatrix, resGraph2$A, obs2$latentAdj)

