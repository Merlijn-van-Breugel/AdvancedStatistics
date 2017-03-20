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
# install.packages("robustbase")
library(robustbase)   
# install.packages("VIM")
library(VIM)   

## Prepare data  
# Choose all relevant columns
cols <- c("Fuel","Price","Cylinders","Displacement","DriveWheel","BHP","Torque","Acceleration",
             "TopSpeed","Weight","Length","Width","Height","AdjustableSteering","Automatic","ClimateControl",
             "CruiseControl","ElectricSeats","ParkingSensors","PowerSteering","SatNav","ESP","Origin")
# We see quite a large number of missing values. Drop all observations for which MPG is NA
TopGear.noNA.MPG <- TopGear[is.na(TopGear$MPG)==FALSE,]
# Choose relevant columns
x <- TopGear.noNA.MPG[,names(TopGear.noNA.MPG) %in% cols]
MPG <- TopGear.noNA.MPG$MPG

# Perform k-Nearest Neighbour imputation method
x.noNA <- kNN(data = x, imp_var = FALSE)


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
    a <- t(apply(z,1, function(x) x/sqrt(sum(x^2))))
    tmp <- matrix(0,ncol(z),ncol(z))
    for(i in 1:nrow(z)){
        tmp <- tmp + a[i,] %*%  t(a[i,])
    }
    return((1/nrow(z)) * tmp)
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
    # D <- diag(apply(z,2,function(x) Qn(x))) 
    # Z <- t(apply(z,1, function(x) solve(D)%*%x))
    # p <- ncol(z)
    # U <- diag(p)
    # for (i in 2:p) {
    #     for (j in 1:(i - 1)) {
    #         U[i, j] <- 1/4*(Qn(Z[,i]+Z[,j])^2 - Qn(Z[,i]-Z[,j])^2) }
    # }
    # E <- eigen(U, symmetric=TRUE)$vectors
    # V <- Z%*%E 
    # L <- diag(apply(V,2,function(x) Qn(x))^2)
    # SIG_rawogk <- D%*%E%*%L%*%t(E)%*%t(D)
    SIG_rawogk <- covOGK(Z,sigmamu= s_Qn)$cov
    return(SIG_rawogk)
}

# sigOGK <- covOGK(Z,sigmamu= s_Qn) #slightly different because of iterations..

# Put all covariance functions in a list
cov.fun.list <- list(corHT = corHT, 
                     corSpearman = corSpearman, 
                     corNSR = corNSR, 
                     covMSS = covMSS, 
                     covBACON1 = covBACON1, 
                     rawCovOGK = rawCovOGK)

# Standardize functions for list
medianWithoutNA<-function(x) {
    median(x[which(!is.na(x))])
}
QnWithoutNA<-function(x) {
    Qn(x[which(!is.na(x))])
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
covDetMCD <- function(x, h0, h, cov.fun.list, eps, delta, ...) {
    
    rownames(x) <- as.character(seq(1:nrow(x)))
    median  <- apply(x, 2, medianWithoutNA)
    Qn      <- apply(x, 2, QnWithoutNA)
    x.dm    <- t(apply(x, 1, function(x) x - median)) 
    Z       <- t(apply(x.dm, 1, function(x) x / Qn)) 
    n.parm  <- NCOL(Z)
    n.obs   <- NROW(Z)    
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
    names(H0)   <- names(cov.fun.list)
    
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
                    det.diff = det.diff,
                    det      = det(H.cov.new),
                    H.cov    = H.cov.new,
                    iter     = iter,
                    H.mean   = colMeans(H),
                    best     = ind))
    }
    
    local.opt <- lapply(H0, function(x) loopTillConv(x,eps))   
    # Get determinants of local.opt
    local.det <- lapply(local.opt, `[[`, 2)
    #print(names(which(local.det==min(unlist(local.det)))))
    opt.cov.fun   <-  names(which(local.det==min(unlist(local.det))))[1]
    # Get rawDetMCD by choosing smallest from local.det and get corresponding output elements from convergence FUN
    rawDetMCD <- local.opt[[opt.cov.fun]]

    #print(rawDetMCD$H.cov)
    #print(rawDetMCD$H.mean)
    
    # Correct the covariance estimate with Fisher correction
    alpha   <- h / nrow(Z)
    cor.fac <- alpha / pgamma(q = qchisq(p = alpha, df = n.parm)/2, shape = (n.parm/2 + 1), rate = 1) 
    # We may want to add a small sample correction
    
    # Get distance with adjusted covariance matrix, using cor.fac
    D       <- mahalanobis(Z, rawDetMCD$H.mean, rawDetMCD$H.cov*cor.fac)
    # Determine weights
    thres   <- qchisq(p = 1-delta, df = n.parm)
    weights <- as.numeric(D <= thres)
    
    # MCD estimates (these are calculated with the non-standardized data x)
    T.rwgt  <- (1/sum(weights)) * t(weights) %*% x
    x.cent.weighted <-  weights * t(apply(x, 1, function(x) x - T.rwgt))
    S.rwgt  <- (1/sum(weights)) * t(x.cent.weighted) %*% x.cent.weighted

    cor.fac <- (1 - delta) / pgamma(q = qchisq(p = 1 - delta, df = n.parm)/2, shape = (n.parm/2 + 1), rate = 1)     
    Sigma.rwgt <- cor.fac * S.rwgt
    
    # Make output list
    output.list <- list(center     = T.rwgt,
                        cov        = Sigma.rwgt,
                        weights    = weights,
                        raw.center = rawDetMCD$mean,
                        raw.cov    = rawDetMCD$cov,
                        best       = rawDetMCD$best,
                        opt.cov.fun = opt.cov.fun)
                        
}

# Test covDetMCD function
test <- covDetMCD(x, h0, h, cov.fun.list, eps, delta) 

## Function for regression based on the deterministic MCD

# Input:
# x .............. matrix of explanatory variables
# y .............. response variable
# MCD.varlist .... containing all setup parameters for covDetMCD
#    h0 .......... size of conservative first subset
#    h ........... subset size
#    cov.fun.list. list with covariance functions
#    eps ......... tolerance level for convergence           
#    delta ....... quantile size for weight determination        

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMCD <- function(x, y, MCD.varlist) {
    MCD   <- covDetMCD(x            = cbind(x,y), 
                       h0           = MCD.varlist$h0, 
                       h            = MCD.varlist$h, 
                       cov.fun.list = MCD.varlist$cov.fun.list, 
                       eps          = MCD.varlist$eps, 
                       delta        = MCD.varlist$delta)
    Sig   <- MCD$cov
    mu    <- MCD$center
    
    # print(Sig)
    # print(mu)
    
    n.parm<- NCOL(x)
    SigXX <- Sig[1:n.parm,1:n.parm]
    SigXY <- Sig[1:n.parm,n.parm+1]
    beta  <- solve(SigXX,SigXY)
    alpha <- drop(MCD$center[n.parm+1] - MCD$center[1:n.parm] %*% beta) # Make scalar instead of matrix
    if (n.parm == 1) { # Can be improved
        yfit  <- alpha + x * beta
    }else{
        yfit  <- alpha + x %*% beta
    }
    res   <- y - yfit
    coefficients <- matrix(c(alpha,beta),1,length(beta)+1)
    colnames(coefficients) <- c("(intercept)",sprintf("beta",seq(1:length(beta))))
    results <- list(coefficients  = coefficients, 
                    fitted.values = yfit, 
                    residuals     = res, 
                    MCD           = MCD )
    return(results)
}

# Multiple Regression function
MultipleRegression <- function(x, y, MCD.varlist, alpha, ...){
    MM   <- lmrob(y~x)
    lm   <- lm(y~x)
    LTS  <- ltsReg(y~x, alpha = alpha)
    MCD  <- lmDetMCD(x,y,MCD.varlist = MCD.varlist) 
    return(list(MCD = MCD,
                lm  = lm,
                LTS = LTS,
                MM  = MM))
    
}
    
# SIMULATION FUNCTION
RobustSimulation <- function(n, sim,  cont.level, cont.type, MCD.varlist, ...){
    coef.list <- list(MCD = matrix(0, sim, 2),
                      lm  = matrix(0, sim, 2),
                      LTS = matrix(0, sim, 2),
                      MM  = matrix(0, sim, 2))
    for(i in 1:sim){
        x <- rnorm(n, mean=5, sd= 1)
        y <- rnorm(n, -2 + x*3, sd=1) 
        if (cont.level > 0){
            if (cont.type == "GoodLev"){
                x[1:(cont.level*n)] <- rnorm(cont.level*n ,mean=15,sd=1)
                y[1:(cont.level*n)] <- rnorm(cont.level*n, -2 + 3*x[1:(cont.level*n)], sd=1)
            }  else if (cont.type == "BadLev"){
                x[1:(cont.level*n)] <- rnorm(cont.level*n ,mean=15,sd=1)
                y[1:(cont.level*n)] <- rnorm(cont.level*n, 5 , sd=1)
            }  else if (cont.type == "VertOut"){
                y[1:(cont.level*n)] <-   rnorm(cont.level*n, 40 , sd=2)
            }
        }
        
        res <- MultipleRegression(x, y, MCD.varlist, alpha = alpha)
        coef.list$MCD[i,]   <- res$MCD$coefficients
        coef.list$lm[i,]    <- res$lm$coefficients
        coef.list$LTS[i,]   <- res$LTS$coefficients
        coef.list$MM[i,]    <- res$MM$coefficients
        
    }
    return(coef.list)
}

# Easier to work with alphas instead of h and h0
n       <- 400
alpha0  <- 0.5
alpha   <- 0.75

MCD.varlist <- list(
                    h0           = alpha0*n,  
                    h            = alpha*n,
                    cov.fun.list = cov.fun.list, 
                    eps          = 1e-10, 
                    delta        = 0.025)

sim1        <- RobustSimulation(n=n, sim=200, cont.level = 0.3, cont.type = "GoodLev", MCD.varlist = MCD.varlist, alpha = alpha)

## Perform grid of simulation specification 
# Contamination type - 3 options ... {Good Leverage, Bad Leverage, Vertical Outlier}
# Contamination level - 6 .......... {0,0.1,...,0.5}
# Subset size - 2 options .......... {0.5,0.75}

cont.types.range <- c("GoodLev","BadLev","VertOut")
cont.level.range <- seq(0,0.5,0.1)
alpha.range      <- c(0.5,0.75)

n                <- 100
sim              <- 1000

MCD.varlist <- list(
    h0           = alpha0*n,  
    h            = NULL,
    cov.fun.list = cov.fun.list, 
    eps          = 1e-10, 
    delta        = 0.025)

SimulationGrid <- function(n, sim,  cont.level.range, cont.types.range, alpha.range, MCD.varlist, ...){
    # Create data frame to write output into
    n.options   <- length(cont.level.range)*length(cont.types.range)*length(alpha.range)
    output      <- as.data.frame(matrix(0,n.options,19)) # 19 comes from 4 estimators, all returning 2 coefficients, mean and var
    counter     <- 1
    for (i in 1:length(alpha.range)){
        MCD.varlist$h <- alpha.range[i]*n
        for (j in 1:length(cont.types.range)){
            for (k in 1:length(cont.level.range)){  
                run <- RobustSimulation(n = n, 
                                        sim = sim, 
                                        cont.level = cont.level.range[k], 
                                        cont.type = cont.types.range[j], 
                                        MCD.varlist = MCD.varlist, 
                                        alpha = alpha.range[i])       
                mean <- lapply(run, function(x) colMeans(x))
                var  <- lapply(run, function(x) sqrt(diag(var(x))))
                if (counter == 1){
                    names(output) <- c("alpha","Cont.Type","Cont.Level",
                                       paste(names(unlist(mean)),"mean", sep = "."),
                                       paste(names(unlist(var)),"var", sep = "."))
                }
                output[counter,(1:3)]  <- cbind(alpha.range[i],cont.types.range[j],cont.level.range[k])
                output[counter,-(1:3)] <- cbind(unlist(mean),unlist(var))
                counter <- counter + 1
            }
        }
    }
    return(output)
}

sim.res <- SimulationGrid(n = n, 
                       sim = sim,  
                       cont.level.range = cont.level.range, 
                       cont.types.range = cont.types.range, 
                       alpha.range= alpha.range, 
                       MCD.varlist = MCD.varlist)

## Do a variable selection on x.noNA
# First make model matrix
x.exp <- model.matrix(~., data = x.noNA)[,-1]

variableSelect <- function(data){
    all.sign <- FALSE
    while (all.sign == FALSE){
        MM.res   <- lmrob(MPG~data, k.max = 500, refine.tol = 1e-05)
        p.val    <- summary(MM.res)$coefficients[,4] 
        if (max(p.val) <= 0.05){
            all.sign <- TRUE
            # summary(MM.res)
        }else{
            drop.var <- gsub("data", "", names(sort(p.val, decreasing = TRUE)[1]))
            var.sel  <- colnames(data)[!(colnames(data) %in% drop.var)] 
            data    <- data[,colnames(data) %in% var.sel]
        }
    }
    return(list(var.sel  = var.sel,
                data.sel = data))
}

data.select <- variableSelect(x.exp)

# Apply multiple regression estimators on dataset with selectedd variables
MCD.varlist <- list(
    h0           = alpha0*nrow(x.exp),  
    h            = alpha*nrow(x.exp),
    cov.fun.list = cov.fun.list, 
    eps          = 1e-10, 
    delta        = 0.025)
TopGear.res <- MultipleRegression(x = data.select$data.sel, y = MPG, MCD.varlist, alpha)

lapply(cov.fun.list, function(x) x(scale(test)))

View(data.select$data.sel)

