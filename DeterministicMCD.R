# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

## Use this code skeleton to implement the deterministic MCD and the resulting
## robust regression estimator.  Please submit your implementation together
## with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will allow me to do some automated testing of your
##            implementations.  Of course you can define any additional
##            functions you may need.

## Load packages
install.packages("robustbase")
library(robustbase)


## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
    a <- apply(z,2, function(x) tanh(x),na.omit=TRUE)
    return(cor(a))
}

# spearman correlation matrix
corSpearman <- function(z) {
    a <- apply(z,2,function(x) rank(x))
    return(cor(a))
}

# correlation matrix based on normal scores of the ranks
corNSR <- function(z) {
    a <- apply(z,2, function(x) qnorm((rank(x)-1/3)/(length(x)+1/3),0,1) )
    return(cor(a))
}

# modified spatial sign covariance matrix
covMSS <- function(z) {
  # *enter your code here*
}

# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  # *enter your code here*
}

# raw OGK estimator of the covariance matrix with median and Qn
# Hint: have a look at function covOGK() in package robustbase
rawCovOGK <- function(z) {
  # *enter your code here*
}



## Main function for deterministic MCD algorithm

# Input:
# x ... data matrix
# h ... subset size
# anything else you need

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# any other output you want to return
covDetMCD <- function(x, h, ...) {
  
    # Standardize data
    medianWithoutNA<-function(x) {
        median(x[which(!is.na(x))])
    }
    QnWithoutNA<-function(x) {
        Qn(x[which(!is.na(x))])
    }
    median <- apply(x, 2, medianWithoutNA)
    Qn <- apply(x, 2, QnWithoutNA)
    x.dm <- t(apply(x, 1, function(x) x - median)) 
    Z <- t(apply(x.dm, 1, function(x) x / Qn)) 
    
    # Get starting values from 6 robust covariance estimates of Z
    # Make list of functions, such that we can elegantly compute Sk
    # Step 1
    covf.list <- list(corHT,corSpearman,corNSR,covMSS,covBACON1,rawCovOGK)
    S.list <- list()
    S.list$S1 <- corHT(Z)
    S.list$S2 <- corSpearman(Z)
    S.list$S3 <- corNSR(Z)
    S.list$S4 <- covMSS(Z)
    S.list$S5 <- covBACON1(Z)
    S.list$S6 <- rawCovOGK(Z)
    
    eigen       <- lapply(S.list, eigen) 
    eigenvec    <- lapply(eigen, `[[`, 2)
    V           <- lapply(eigenvec, function(x) Z %*% x) 
    
    # Step 2
    Qn.V        <- lapply(V, function(x) apply(x, 2, QnWithoutNA))
    L           <- lapply(Qn.V, function(x) diag(x^2))
    Sigma       <- lapply(seq(1,length(S.list)), function(x) eigenvec[[x]] %*% L[[x]] %*% t(eigenvec[[x]]))
    
    # Step 3
    Sigma.chol  <- lapply(Sigma,function(x) chol(x))
    mu          <- lapply(Sigma.chol,function(x) x %*% apply(Z %*% solve(x), 2, medianWithoutNA))
    
    # Get distance 
    D           <- lapply(seq(1,length(S.list)), function(x) mahalanobis(Z, mu[[x]],Sigma[[x]]))
}



## Function for regression based on the deterministic MCD

# Input:
# x ... matrix of explanatory variables
# y ... response variable
# h ... subset size
# anything else you need

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMDC <- function(x, y, h, ...) {
  # *enter your code here*
}
