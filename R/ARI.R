#' Evalute the adjusted Rand index
#'
#' Compute the adjusted Rand index to compare two partitions
#'
#' @details
#' the partitions may be provided as n-vectors containing the cluster memeberships
#'     of n entities, or by Qxn - matrices whose entries are all
#'     0 and 1 where 1 indicates the cluster membership
#' @param x vector (of length n) or matrix (with n columns) providing a partition
#' @param y vector or matrix providing a partition
#'
#' @return the value of the adjusted Rand index
#' @export
#'
#' @examples
#' clust1 <- c(1,2,1,2)
#' clust2 <- c(2,1,2,1)
#' ARI(clust1, clust2)
#'
#' clust3 <- matrix(c(1,1,0,0, 0,0,1,1), nrow=2, byrow=TRUE)
#' clust4 <- matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,1), nrow=3, byrow=TRUE)
#' ARI(clust3, clust4)
ARI <- function(x, y) {
  if (is.matrix(x))
    x <- apply(x,2,which.max)
  if (is.matrix(y))
    y <- apply(y,2,which.max)
  # first, get crosstabs
  ctab <- table(x,y)
  # now calculate 4 intermediary sums
  cellsum <- sum(ctab*(ctab-1)/2)
  totsum <- sum(ctab)*(sum(ctab)-1)/2
  # use matrix multiplication to get row and column marginal sums
  rows <- ctab %*% rep(1,ncol(ctab))
  rowsum <- sum(rows*(rows-1)/2)
  cols <- rep(1,nrow(ctab)) %*% ctab
  colsum <- sum(cols*(cols-1)/2)
  # now put them together
  adj.rand <- (cellsum - (rowsum*colsum/totsum))/(.5*(rowsum +colsum)-(rowsum*colsum/totsum))
  return(adj.rand)
}
