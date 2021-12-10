
#' @title Computes the non-parametric estimator of the multivariate Spearman's rho measure of average upper orthant dependence proposed by Pérez and Prieto-Alaiz (2016).
#' @description Computes the non-parametric estimator of the multivariate Spearman's rho measure of average upper orthant dependence proposed by Pérez and Prieto-Alaiz (2016).
#' @param X nxd matrix of the data
#' @return value of the coefficient
#' @export rhoplus_PP
#' @examples
#' library(copula)
#' set.seed(111)
#' mv.Cl.NN <- mvdc(claytonCopula(dim = 3, iRho(claytonCopula(),0.7)), c("norm", "norm", "norm"),
#' list(list(mean = 0, sd =2), list(mean = 0, sd =2), list(mean = 0, sd =2)))
#' x.sample <- rMvdc(10000, mv.Cl.NN)
#' rhoplus_PP(x.sample)
rhoplus_PP <- function (X)
{
  ## X is the nxd matrix of our data
  ## n=sample size
  ## d=number of dimensions
  library(copula)
  R <- pobs(X)*(length(X[,1])+1) #the nxd matrix of ranks
  d<-ncol(R)
  n<-nrow(R)
  c0<-mean(apply(R,1,prod))
  c1<-((n+1)/2)^d
  c2<-mean((seq(1:n))^d)
  value <- round((c0-c1)/(c2-c1),3)
  return(value)
}


#' @title Computes the non-parametric estimator of the multivariate Spearman's rho measure of average lower orthant dependence proposed by Pérez and Prieto-Alaiz (2016).
#' @description Computes the non-parametric estimator of the multivariate Spearman's rho measure of average lower orthant dependence proposed by Pérez and Prieto-Alaiz (2016).
#' @param X nxd matrix of the data
#' @return value of the coefficient
#' @export rhominus_PP
#' @examples
#' library(copula)
#' set.seed(111)
#' mv.Cl.NN <- mvdc(claytonCopula(dim = 3, iRho(claytonCopula(),0.7)), c("norm", "norm", "norm"),
#' list(list(mean = 0, sd =2), list(mean = 0, sd =2), list(mean = 0, sd =2)))
#' x.sample <- rMvdc(10000, mv.Cl.NN)
#' rhominus_PP(x.sample)
rhominus_PP <- function (X)
{
  ## X is the nxd matrix of our data
  ## n=sample size
  ## d=number of dimensions
  library(copula)
  R <- pobs(X)*(length(X[,1])+1) #the nxd matrix of ranks
  d<-ncol(R)
  n<-nrow(R)
  Rbar<-(n+1)-R
  c0<-mean(apply(Rbar,1,prod))
  c1<-((n+1)/2)^d
  c2<-mean((seq(1:n))^d)
  value <- round((c0-c1)/(c2-c1),3)
  return(value)
}

#' @title Computes the tie-corrected non-parametric estimator of the multivariate Spearman's rho measure of average lower orthant dependence proposed by Genest et al. (2013).
#' @description Computes the tie-corrected non-parametric estimator of the multivariate Spearman's rho measure of average lower orthant dependence proposed by Genest et al. (2013).
#' @param data nxd matrix of the data
#' @return value of the coefficient
#' @export rhominus_GNR
#' @examples
#' library(copula)
#' set.seed(111)
#' mv.Cl.NN <- mvdc(claytonCopula(dim = 3, iRho(claytonCopula(),0.7)), c("norm", "norm", "norm"),
#' list(list(mean = 0, sd =2), list(mean = 0, sd =2), list(mean = 0, sd =2)))
#' x.sample <- rMvdc(10000, mv.Cl.NN)
#' rhominus_GNR(x.sample)
rhominus_GNR <- function(data){ # data= nxd matrix of data
  X <- data
  library(copula)
  R <- pobs(X)*(length(X[,1])+1) #the nxd matrix of ranks
  #X<-X[,1:3]
  d <- ncol(R)
  n <- nrow(R)
  constant <- (d + 1)/((2^d)-d-1)
  value <- round(-constant + constant*(2^d)*mean(apply(((2*n+1)/(2*n))-(R/n), 1, prod)), 3)
  return(value)
}


#' @title Computes the tie-corrected non-parametric estimator of the multivariate Spearman's rho measure of average upper orthant dependence proposed by Genest et al. (2013).
#' @description Computes the tie-corrected non-parametric estimator of the multivariate Spearman's rho measure of average upper orthant dependence proposed by Genest et al. (2013).
#' @param data nxd matrix of the data
#' @return value of the coefficient
#' @export rhoplus_GNR
#' @examples
#' library(copula)
#' set.seed(111)
#' mv.Cl.NN <- mvdc(claytonCopula(dim = 3, iRho(claytonCopula(),0.7)), c("norm", "norm", "norm"),
#' list(list(mean = 0, sd =2), list(mean = 0, sd =2), list(mean = 0, sd =2)))
#' x.sample <- rMvdc(10000, mv.Cl.NN)
#' rhoplus_GNR(x.sample)
rhoplus_GNR <- function(data){ #data = nxd matrix of data
  X <- data
  library(copula)
  R <- pobs(X)*(length(X[,1])+1) #the nxd matrix of ranks
  #X <- X[,1:3]
  d <- ncol(R)
  n <- nrow(R)
  constant <- (d + 1)/((2^d)-d-1)
  value <- round(-constant + constant*(2^d)*mean(apply(((R/n) - (1/(2*n))), 1, prod)), 3)
  return(value)
}

