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

## Prepare data

# We see quite a large number of cells with NA value
# This is problematic for correlation and distance calculation, hence we omit these rows
x <- na.omit(TopGear[,9:15])

## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
    a <- apply(z,2, function(x) tanh(x))
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
    a <- apply(z,1, function(x) x/sqrt(sum(x^2)))
    return((1/nrow(z))*a%*%t(a))
}

# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
    a <- apply(z,1,function(x) sqrt(sum(x^2)))
    n <- ceiling(nrow(z)/2)
    b <- z[order(a)[1:n],]
    return(cov(b))
}

# raw OGK estimator of the covariance matrix with median and Qn
# Hint: have a look at function covOGK() in package robustbase
rawCovOGK <- function(z) {
    D <- diag(apply(z,2,function(x) Qn(x))) 
    Z <- apply(z,1, function(x) solve(D)%*%x)
    p <- ncol(z)
    U <- diag(p)
    for (i in 2:p) {
        for (j in 1:(i - 1)) {
            U[i, j] <- 1/4*(Qn(Z[,i]+Z[,j])^2 - Qn(Z[,i]-Z[,j])^2)}
    }
    E <- eigen(U, symmetric=TRUE)$vectors
    V <- Z%*%E 
    L <- diag(apply(V,2,function(x) Qn(x))^2)
    SIG_rawogk <- D%*%E%*%L%*%t(E)%*%t(D)
    return(SIG_rawogk)
}

# Put all covariance functions in a list
cov.fun.list <- list(corHT, corSpearman, corNSR, covMSS, covBACON1, rawCovOGK)
cov.fun.list <- list(corHT, corSpearman, corNSR) 

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
covDetMCD <- function(x, h0, h, cov.fun.list, eps, delta, ...) {
    
    # Standardize data
    medianWithoutNA<-function(x) {
        median(x[which(!is.na(x))])
    }
    QnWithoutNA<-function(x) {
        Qn(x[which(!is.na(x))])
    }
    median  <- apply(x, 2, medianWithoutNA)
    Qn      <- apply(x, 2, QnWithoutNA)
    x.dm    <- t(apply(x, 1, function(x) x - median)) 
    Z       <- t(apply(x.dm, 1, function(x) x / Qn)) 
    n.parm  <- ncol(Z)
    n.obs   <- nrow(Z)    
    # Get starting values from 6 robust covariance estimates of Z
    # Make list of functions, such that we can elegantly compute Sk
    # Step 1
    S.list <- lapply(cov.fun.list, function(x) x(Z)) 
    
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
   
    # Get h smallest distance values
    ind         <- lapply(D, function(x) sort(x)[1:h0]) 
    H0          <- lapply(ind, function(x) Z[names(x),]) 
    
    
    ## Write function to converge to smallest subset
    loopTillConv <- function(H, eps){        
        det.diff <- Inf
        iter <- 1
        while (det.diff > eps){
            
            # Compute mean and covariance of subset of Z (H)
            H.mean  <- colMeans(H)
            H.cov   <- cov(H)
            # Get new mahalanobis distance of all data points Z, based on H.mean and H.cov
            D       <- mahalanobis(Z, H.mean,H.cov)
            # Get h points with smallest distance
            ind     <- sort(D)[1:h]
            H       <- Z[names(ind),] 
            H.cov.new <- cov(H)
            
            det.diff <- det(H.cov) - det(H.cov.new)
            iter    <- iter + 1
        }
        return(list(
                    det.diff <- det.diff,
                    det   <- det(H.cov.new),
                    H.cov <- H.cov.new,
                    iter  <- iter,
                    H.mean <- colMeans(H)),
                    best  <- ind)
    }
    
    local.opt <- lapply(H0, function(x) loopTillConv(x,eps))   
    # Get determinants of local.opt
    local.det <- lapply(local.opt, `[[`, 2)
    # Get rawDetMCD by choosing smallest from local.det and get corresponding output elements from convergence FUN
    rawDetMCD <- local.opt[which(local.det==min(unlist(local.det)))[1]]
    
    # Correct the covariance estimate with Fisher correction
    alpha   <- h / nrow(Z)
    cor.fac <- alpha / pgamma(q = qchisq(p = alpha, df = n.parm)/2, shape = (n.parm/2 + 1), rate = 1) 
    # We may want to add a small sample correction
    
    # Get distance with adjusted covariance matrix, using cor.fac
    D       <- mahalanobis(Z, rawDetMCD[[1]], rawDetMCD[[3]]*cor.fac)
    # Determine weights
    thres   <- qchisq(p = 1-delta, df = n.parm)
    weights <- as.numeric(D <= thres)
    
    # MCD estimates
    T.rwgt  <- (1/sum(weights)) * t(weights) %*% Z
    Z.cent.weighted <-  weights * t(apply(Z, 1, function(x) x - T.rwgt))
    S.rwgt  <- (1/sum(weights)) * t(Z.cent.weighted) %*% Z.cent.weighted

    cor.fac <- (1 - delta) / pgamma(q = qchisq(p = 1 - delta, df = n.parm)/2, shape = (n.parm/2 + 1), rate = 1)     
    Sigma.rwgt <- cor.fac * S.rwgt
    
    # Make output list
    output.list <- list(center     <- T.rwgt,
                        cov        <- Sigma.rwgt,
                        weights    <- weights,
                        raw.center <- rawDetMCD$mean,
                        raw.cov    <- rawDetMCD$cov,
                        best       <- rawDetMCD$best)
                        
                        
    )
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
