# --------------------------------------------------------------------
# Author: Merlijn van Breugel
# Date: 18 April 2017
# --------------------------------------------------------------------

################################################################################
## Disclaimer from lecturer                                                    #                 
################################################################################
## Use this code skeleton to implement the procedures for obtaining point      #
## estimates and valid standard errors via multiple imputation and the         #
## nonparametric bootstrap.  Please submit your implementation together        #
## with your report.                                                           #
##                                                                             #
## IMPORTANT: Please do not change any function names and make sure that you   #
##            use the correct input and output as specified for each function. #
##            This will simplify grading because each student's code will be   #
##            similarly structured.  However, feel free to add other arguments #
##            to the function definitions as needed.                           #    
################################################################################

################################################################################ 
## Initialization                                                              #
################################################################################

## Setup of R
# Set working directory
setwd("C:\\Users\\Merlijn\\Documents\\GitHub\\AdvancedStatistics")

# Load packages
require(foreach)
require(parallel)
require(doParallel)
require(snow)
require(VIM)
require(robustbase)  
require(MASS)

# Set seed. Best-practice: large seed for more randomness
set.seed(42424242)

# Make sure functions are loaded, needed for parallel function
source('createBaseData.R')
source('createSimData.R')
source('multimp.R')
source('fit.R')
source('pool.R')


################################################################################
## Create functions                                                            #                 
################################################################################

## Functions for multiple imputation via iterative model-based imputation

# Multiple imputation
# Input:
# x ............ data set with missing values
# m ............ number of imputations
# Additional arguments to be passed to function irmi() from package VIM
# robust ....... whether to use robust methods
# robMethod .... which robust method to use for continuous response ...
#                if 'None' then no Robust Method is applied
# Output:
# A list containing the m imputed data sets
multimp <- function(x, m = NULL, mixed = NULL, robust = FALSE, robMethod = 'MM', 
                    init.method = 'kNN',...) {
    
    # Set a sensible default for the number of imputations m
    # Follow advice of Bodner (2008) and White et al. (2011), set m default to.. 
    # ... percentage of missing information
    if (is.null(m)){
        m <- 5*round(100*mean(is.na(x))/5)
    }
    ## Use function irmi() from VIM package multiple imputations
    xList <- irmi(x, 
                  maxit = 100, 
                  mi = m, 
                  mixed = mixed,
                  init.method = init.method,
                  robust = robust,
                  robMethod = robMethod) 
    names(xList) <- paste0('Imputation',1:m)
    return(xList) 
}

# Fit regression models
# Input:
# xList .. list of imputed data sets as returned by function multimp()
# y.name . name of column to use as y variable 
# method . character string specifying whether to use OLS or MM-estimator
# ........ additional arguments to be passed to modeling function (for example,
# ........ control parameters for the MM-algorithm)
# Output:
# A list of regression models (each fitted to one of the imputed data sets)
fit <- function(xList, y.name = 'y', method = 'MM', method.inputList = NA, ...){
    
    ## Perform selected regression method
    # Create formula for regression
    fm <- paste0(y.name,'~',paste0(names(xList[[1]])[(names(xList[[1]]) != y.name)],collapse = '+'))
    # Fit regression models
    if (method == 'MM'){
        fit <- lapply(xList, function(x){lmrob(fm,data = x, method = 'MM',
                                               setting = 'KS2014',
                                               k.max = 300)})
    }else if (method == 'OLS'){
        fit <- lapply(xList, function(x){lm(fm,data = x)})
    }
    
    return(fit) 
}


# Pool point estimates and standard errors
# Input:
# fitList  list of fitted regression models as returned by function fit()
# Output:
# A matrix that contains the pooled coefficients in the first column, the
# pooled standard errors in the second column, the t-statistic in the third
# column, the estimated degrees of freedom in the fourth column, and the
# p-value in the fifth column (see slide 50 of Lecture 5 for an example)
pool <- function(fitList, return.more = FALSE, method = 'OLS',...){
    
    # Retrieve useful parameters
    names.var  <- names(fitList[[1]][['coefficients']])
    n.obs <- NROW(fitList[[1]]$x)
    n.var  <- length(fitList[[1]]$coefficients)
    m <- length(fitList)

    # Get converged imputations if MM is used, OLS works fine
    if (method == 'MM'){
        fitList <- lapply(fitList,function(x){
                      if (x[['converged']] == TRUE){x}})
        fitList <- fitList[!sapply(fitList, is.null)] 
    }
    
    # Start with pooled coefficients estimates
    mean.pool <- sapply(seq(1:n.var), function(i)
    {
        coef <- sapply(fitList, function(x){coef(summary(x))[i,'Estimate']})
        mean(coef)
    }) 
    
    # Compute within and between inputed variance
    U.var <- sapply(seq(1:n.var), function(i)
    {
        sd <-   sapply(fitList, function(x){coef(summary(x))[i,'Std. Error']})
        mean(sd * sd)
    })
    B.var <- sapply(seq(1:n.var), function(i)
    {
        b <-   sapply(fitList, 
                      function(x){
                          (coef(summary(x))[i,'Estimate']
                           - mean.pool[i])^2})
        (1/(m-1)) * sum(b)
    })
    # Combine both variances into pooled variance
    var.pool <- U.var + ((m+1)/m) * B.var 
    sd <- sqrt(var.pool)
    
    # Get degrees of freedom, following Barnard and Rubin (1999) 
    df.comp <- n.obs-n.var-1 # QUESTION: Do you count the intercept as explanatory variable? 
    frac.miss <- ((m+1)/m)*(B.var/var.pool)
    df.obs <- ((df.comp+1)/(df.comp+1))*df.comp*(1-frac.miss)
    df.m <- (m-1)/(frac.miss^2)
    df <- (df.m*df.obs)/(df.m+df.obs)
    
    # Get t-statistics based on pooled mean and variance
    t.stat <- mean.pool/sd
    # Two-sided, get p value
    p.val <- 2*pt(abs(t.stat), df = df, lower = FALSE) 
    # Nice function to get significance
    symp <- symnum(p.val, corr = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))
    
    output.table <- (data.frame('Estimate' = mean.pool,
                                         'Std.Error' = sd,
                                         't.value' = t.stat,
                                         'df' = df,
                                         'p.val' = p.val,
                                         'sign' = symp ))
    rownames(output.table) <- names.var    
    # If asked by user, also return average robustness weights
    if (return.more == TRUE){
        mean.rweights <- sapply(seq(1:n.obs), function(i)
        {
            rweights <- sapply(fitList, function(x){x[['rweights']][i]})
            mean(rweights)
        })        
        output <- list(output.table = output.table,
                       mean.rweights = mean.rweights,
                       imp.sets = length(fitList))
    }else{
        output <- output.table
    }
    
    return(output)
}

## Function for the bootstrap with kNN imputation and linear regression
# Input:
# x       data set with missing values
# R       number of replications
# k       number of neighbors for kNN imputation
# method  character string specifying whether to use OLS or MM-estimator
# Output:
# A list with the following components:
# replicates  a matrix with all coefficient estimates from all replications
# summary     a matrix that contains the point estimates of the coefficients in
#             the first column, the standard errors in the second column, the
#             z-statistic in the third column, and the p-value in the fourth
#             column (see slide 29 of Lecture 5 for an example)
bootstrap <- function(x, y.name = 'y', R = 1000, k = 5, method = 'MM', ...) {
    
    # Retrieve useful parameters
    names.var  <- names(x)
    n.obs <- NROW(x)
    n.var <- (NCOL(x)-1)
    # Already construct formula for regression
    fm <- paste0(y.name,'~',paste0(names(x)[(names(x) != y.name)],collapse = '+'))
    
    ## Get R bootstrap samples with regression results
    # Massive improvement by parallelization
    cores <- detectCores()
    cl <- makePSOCKcluster(cores)
    registerDoParallel(cl)
    # Ugly, but needed: all workers need to have the packages and functions to their disposal
    clusterCall(cl, function(){library(VIM)})
    clusterCall(cl, function(){library(robustbase)})

    # time.start <- Sys.time()
    bootstrap.res <-  foreach(i=1:R, .combine=rbind) %dopar% {  
        
        s <- x[sample(NROW(x), n.obs, replace = TRUE),]
        # We use kNN method from the VIM package for single imputation
        sample.imp <- kNN(s, k = k, imp_var = FALSE)
        # Fit regression models
        if (method == 'MM'){
            fit <- tryCatch( 
                lmrob(fm,data = sample.imp, method = 'MM', setting="KS2014",k.max = 200)
                ,error=function(e){e})
            # we have an instance that is corrupted
            if (any(class(fit) %in% "error")) { 
                R <- R+1 # Perform an extra iteration
            }
        }else if (method == 'OLS'){
            fit <- lm(fm,data = sample.imp)
        }
        fit[['coefficients']] 
    }
    # time.end <- Sys.time()
    # time.end - time.start 
    # Stop parallel cluster workers
    stopCluster(cl)
    coef.mean <- colMeans(bootstrap.res)
    
    # Compute standard deviation of bootstrap replicates
    sd <- sqrt((1/(R-1))*colSums(sweep(bootstrap.res,2,coef.mean)^2))  
    
    # Assume asymptotic normality
    z.stat <- coef.mean/sd

    # Two-sided, get p value
    p.val <- 2*pnorm(abs(z.stat), lower = FALSE) 
    # Nice function to get significance
    symp <- symnum(p.val, corr = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))
    output.table <- (data.frame('Estimate' = coef.mean,
                                         'Std.Error' = sd,
                                         'z.value' = z.stat,
                                         'p.val' = p.val,
                                         'sign' = symp))
    return(output.table)    
}

## Functions to create dataset with missing data and outliers

## Create base data
## Input:
## Data itself: 
# n.obs .............. number of observations
# n.var .............. number of variables, dependent variable y excluded 
# mean.x ............. give means of explanatory variables (for rnorm) 
# factor.cor ......... defines how x[,nvar] is related to x[,nvar-1]
# sd.x ............... give standard deviation of explanatory variables (for rnorm)
# sd.y ............... give standard deviation of dependent variable (for rnorm)
# c .................. intercept for dependent variable
# coef ............... coefficients for dependencies between explanatory variables and dependent variables
## Output:
# A single data set for simulation purposes

createBaseData <- function(n.obs = 100, n.var = 3, 
                           mean.x = c(1,2,3),
                           corx1x2 = 0.5, 
                           sd.x = c(1,1,1), sd.y = 1,
                           c = -5, coef = c(1,2,3)){
    
    ## Create base data set
    # Create p-1 explanatory variables that are normally distributed 
    # We let two variables be correlated to mimick the biopics set 
    
    # Make correlation matrix for explanatories, correlate between x1 and x2
    if (n.var >= 2 & !is.na(corx1x2)){    
    cor <- as.matrix(diag(sd.x),n.var,n.var) 
    cor[1,2] <- corx1x2 * sd.x[1] * sd.x[2]
    cor[2,1] <- corx1x2 * sd.x[1] * sd.x[2]
    x <- mvrnorm(n = n.obs, mu = mean.x, Sigma = cor)
    }else{
        y <- rnorm(n.obs, mean.x[1], sd = sd.x[1])
    }
    
    colnames(x) <- paste0('x',1:n.var)
    # Compute y, given the pre-defined coefficients
    y <- rnorm(n.obs, mean = c + x %*% coef, sd = sd.y)
    
    return(cbind(y,x))
}

## Include missing data and/or outliers
## Input:
## Modification for missings and outliers:
# mechanism.miss...... type of missing data (MCAR, MAR, MNAR,none)
# prop.miss .......... proportion of the data that is missing
# mechanism.outlier... type of outliers within regression context (GoodLev,BadLev,VertOut,None)
# prop.outlier........ proportion of observations that are outliers 
## Output:
# A single data set for simulation purposes

createSimData <- function(base.data, 
                          n.obs = 100, n.var = 3, 
                          mean.x = c(5,10,15), 
                          sd.x = c(1,2,3), sd.y = 2,
                          c = -5, coef,
                          mechanism.miss = 'MCAR', prop.miss = 0.1, 
                          mechanism.outlier = 'BadLev', prop.outlier = 0.1){
    
    y <- base.data[,1]
    x <- base.data[,-1]
    size <- n.obs * (n.var + 1)
    ## Include outliers, re-use part of code from Group Assignment
    # Assign outliers randomly of different explanatory variables
    # Count number of outliers, till prop.outlier achieved
    if (prop.outlier > 0 & mechanism.outlier != 'None'){
        n.outliers <- 0
        outlying.YN <- matrix(0,n.obs,n.var)
        while (n.outliers < prop.outlier * size){
            outlying.var <- sample(1:n.var,1)
            outlying.obs <- sample(1:n.obs,1)
            # Given the outlier type, include outliers
            if (mechanism.outlier == 'GoodLev'){
                x[outlying.obs,outlying.var] <- rnorm(1,
                                                      mean = 0.3*mean.x[outlying.var], 
                                                      sd = 1)
                y[outlying.obs]             <- rnorm(1, 
                                                     mean = c + x[outlying.obs] %*% coef, 
                                                     sd = 1)
            }  else if (mechanism.outlier == 'BadLev'){
                x[outlying.obs,outlying.var] <- rnorm(1,
                                                      mean = 0.3*mean.x[outlying.var], 
                                                      sd = 1)
                y[outlying.obs]             <- rnorm(1, 
                                                     mean = 0.5*mean(y), 
                                                     sd=1)
            }  else if (mechanism.outlier == 'VertOut'){
                y[outlying.obs]             <- rnorm(1, 
                                                     mean = 0.1*mean(y), 
                                                     sd = 1)
            }
            outlying.YN[outlying.obs,outlying.var] <- 1
            n.outliers <- sum(outlying.YN)
        }
    }
    data <- cbind(y,x)
    miss <- matrix(1,n.obs,(n.var+1)) # Introduce missings by NA
    if (!(prop.miss == 0 | mechanism.miss == 'None')){
        ## Replace prop.miss by NA
        if (mechanism.miss == 'MCAR'){
            # Draw random samples from [0,1] uniform distribution
            mcar <- matrix(runif(2*n.obs, min=0, max=1),n.obs,2)
            # All values smaller than the proportion of missings are replaced by NA
            # Introduce missings in both y and x
            miss[,1:2] <- ifelse(mcar>((n.var/2)*prop.miss),1,NA)
            data <- data * miss
        }else if(mechanism.miss == 'MAR'){
            # Create MAR as follows: a covariate is missing, if another variable is larger than its miss.prop quantile
            quantiles <- apply(data,2,function(x){
                quantile(x, seq(from = .01, to = .99, by = .01))})
            # We impute missings in two variables
            # In both roughly equally many missings
            # Hence, the corresponding quantile should be greater than 100-floor((n.var/2)*prop.miss*100)
            miss[,1] <- ifelse(data[,2] < quantiles[100-floor((n.var/2)*prop.miss*100),2],1,NA)                    
            miss[,2] <- ifelse(data[,3] < quantiles[100-floor((n.var/2)*prop.miss*100),3],1,NA)
            data <- data * miss
        }else if(mechanism.miss == 'MNAR'){
            # Values are missing conditional on their own value
            quantiles <- apply(data,2,function(x){
                quantile(x, seq(from = .01, to = .99, by = .01))})
            miss[,1] <- ifelse(data[,1] < quantiles[100-floor((n.var/2)*prop.miss*100),1],1,NA)                    
            miss[,2] <- ifelse(data[,2] < quantiles[100-floor((n.var/2)*prop.miss*100),2],1,NA)
            data <- data * miss
        }
    }
    return(list(
        data = data,
        prop.miss = sum(is.na(data))/size))
}

## Function to perform multiple simulations oevr a variety of configurations
SimulationGrid <- function(sim = 100, base.data.config, 
                           m = NULL,
                           k = 5,
                           mechanism.outlier.range, prop.outlier.range,
                           mechanism.miss.range, prop.miss.range, ...){
    
    # Retrieve some useful parameters
    n.var <- base.data.config$n.var
    n.var.c <- n.var + 1 # Including intercept, easy for notation
    n.obs <- base.data.config$n.obs
    
    # Create data frame to write output into
    n.options   <- length(mechanism.outlier.range)*
        length(mechanism.miss.range)*
        length(prop.outlier.range)*
        length(prop.miss.range)
    n.est  <- n.var.c * 2 * 2 #Single and multiple imputation, with OLS and MM
    res <- array(NA,dim = c(sim,n.est,n.options))
    
    # Start by creating the simulation base data set
    # This split in data generation process is relevant, as we do not want to sample new data for all configuration
    # Given a base set, we want inclusion of NA and outliers
    # We repeat this process for a variety of random base sets, to check for consistency
    
    # Count in which configuration the simulation is, for storing results
    count.config <- 1
    config.name <- matrix(NA,n.options,1)
    
    # Massive improvement by parallelization
    cores <- detectCores()
    cl <- makePSOCKcluster(cores)
    registerDoParallel(cl)
    # Ugly, but needed: all workers need to have the packages and functions to their disposal
    clusterCall(cl, function(){library(VIM)})
    clusterCall(cl, function(){library(robustbase)})
    clusterCall(cl, function(){library(MASS)})
    clusterCall(cl, function(){source('createBaseData.R')})
    clusterCall(cl, function(){source('createSimData.R') })
    clusterCall(cl, function(){source('multimp.R') })
    clusterCall(cl, function(){source('fit.R') })
    clusterCall(cl, function(){source('pool.R') })
    
    # Now loop over all configurations
    for (i in 1:length(mechanism.miss.range)){
        for (j in 1:length(prop.miss.range)){
            for (k in 1:length(mechanism.outlier.range)){
                for (l in 1:length(prop.outlier.range)){
                    # For each configuration, we run this sim times
                    config.run <- foreach(s=1:sim, .combine=rbind) %dopar% {
                        base.data <- createBaseData(n.obs = n.obs, 
                                                    n.var = n.var, 
                                                    mean.x = base.data.config$mean.x, 
                                                    sd.x = base.data.config$sd.x, 
                                                    sd.y = base.data.config$sd.y,
                                                    c = base.data.config$c, 
                                                    coef = base.data.config$coef)
                        
                        sim.data <- createSimData(base.data = base.data,
                                                  n.obs = n.obs,
                                                  n.var = n.var,
                                                  mean.x = base.data.config$mean.x,
                                                  sd.x = base.data.config$sd.x,
                                                  sd.y = base.data.config$sd.y,
                                                  c = base.data.config$c,
                                                  coef = base.data.config$coef,
                                                  mechanism.miss = mechanism.miss.range[i],
                                                  prop.miss = prop.miss.range[j],
                                                  mechanism.outlier = mechanism.outlier.range[k],
                                                  prop.outlier = prop.outlier.range[l])$data
                        coef.names <- c('c',colnames(base.data[,-1]))
                        colnames(sim.data) <- coef.names # Weird, should not be needed
                        regress.coef <- NULL
                        
                        sim.data.kNN <- as.matrix(kNN(sim.data, k = k, imp_var = FALSE))
                        xList <- multimp(x=sim.data, m = NULL)

                        for (method.sim in c('OLS','MM')){
                            # Run single imputation kNN and multiple
                            for (method.imp in c('single','multiple')){
                                if (method.imp == 'single'){
                                    if (method.sim == 'OLS'){
                                        fit.KNN <- lm(sim.data.kNN[,1]~sim.data.kNN[,-1])
                                    }else if (method.sim == 'MM'){
                                        fit.KNN <- lmrob(sim.data.kNN[,1]~sim.data.kNN[,-1], 
                                                         method = 'MM',
                                                         setting = 'KS2014',
                                                         k.max = 200)
                                    }
                                    coef <- fit.KNN$coefficients
                                }else if (method.imp == 'multiple'){
                                    fitList <- fit(xList = xList, y.name = 'y',method = method.sim)
                                    coef <- pool(fitList)[,1]
                                }
                                names(coef) <- paste(method.sim,method.imp,coef.names,sep = '.')
                                regress.coef <- cbind(regress.coef,t(as.matrix(coef)))
                            }
                            
                        }
                        regress.coef
                    }
                    # Store rbinded results into array
                    res[,,count.config] <- config.run
                    config.name[count.config] <- paste0(
                        'MissingMechanism:',mechanism.miss.range[i], 
                        '|PropMissing:',prop.miss.range[j], 
                        '|MechanismOutlier:',mechanism.outlier.range[k], 
                        '|PropOutliers:',prop.outlier.range[l])
                    noquote(paste('Configuration:',count.config))
                    count.config <- count.config + 1
                    
                }
            }
        }
    }
    stopCluster(cl)
    # For convenience, set dimnames
    dimnames(res)[[1]] <- paste('sim',1:sim,sep='.')
    dimnames(res)[[2]] <- names(regress.coef[,-(1:(n.var+1))])
    dimnames(res)[[3]] <- config.name
    
    # Compute mean and variance of coefficients over simulations
    coef.means <- apply(res,2,function(x){colMeans(x)})
    coef.var <- apply(res,2,function(x){sqrt(diag(var(x)))})
    
    # Get quantiles of simulation coefficients
    coef.quantiles <- apply(res,2,function(x){quantile(x, c(.05,.95)) })
    
    output.list <- list(res = res,
                        coef.means = coef.means,
                        coef.var = coef.var,
                        coef.quantiles = coef.quantiles)
    return(output.list)        
}

################################################################################
## Perform simulation study                                                    #                 
################################################################################

# Setup settings of simulation
base.data.config <- list(n.obs = 200, n.var = 3, 
                         mean.x = c(3,4,5),
                         cor = 0.5,
                         sd.x = c(1,1,1), sd.y = 1,
                         c = -5, coef = c(1,2,3))

mechanism.outlier.range <- c("GoodLev","BadLev","VertOut")
mechanism.miss.range <- c("MCAR","MAR","MNAR")
prop.miss.range <- seq(0.1,0.3,0.1)
prop.outlier.range <- seq(0,0.2,0.2)

mechanism.outlier.range <- c("None")
mechanism.miss.range <- c("None","MCAR","MAR","MNAR")
prop.miss.range <- seq(0.1,0.3,0.1)
prop.outlier.range <- 0

## Run simulation
# Warning, this may take quite a long time
set.seed(20170423)

## Please not that both data generating functions should be source-able from file,
# such that the parallel workers can reach them 
sim.res <- SimulationGrid(sim = 10, 
                          base.data.config = base.data.config,  
                          mechanism.outlier.range = mechanism.outlier.range, 
                          prop.outlier.range = prop.outlier.range,
                          mechanism.miss.range = mechanism.miss.range, 
                          prop.miss.range = prop.miss.range)


################################################################################
## Analyze the biopics data                                                    #                 
################################################################################

# Load data
load("biopics.RData")
# Convert variables to numeric whenever possible
# apply(biopics,2,function(x) mode(x))
var.num <- c('year_release','box_office','number_of_subjects','person_of_color')
biopics[,var.num] <- apply(biopics[,var.num],2,function(x) as.numeric(x))

# Transform box office using a log transformation to get more readable estimates
biopics$box_office <- log(biopics$box_office)

# Store names of variable of interest
inc.vars <- c('country','year','box.office','n.subject',
              'type','race','color','gender')

# VIM package does not allow for long variable names in plot labels
# therefore, map variable to shorter ones
biopics.NA <- biopics
names(biopics.NA) <- c('title','country','year','box.office','n.subject',
                       'type','race','color','gender')

################################################################################
## Visualize the missings using VIM package                                    #
## Only code for relevant plots are given                                      #
################################################################################

# Which variables contain missings?
var.mis <- colnames(biopics.NA[,colSums(is.na(biopics.NA))>0])
aggr(biopics.NA[,var.mis], numbers = TRUE,only.miss = TRUE, bars = TRUE)

# Year seems interesting
# should be aggregated for a nice plot
spineMiss(biopics.NA[, c('year','box.office')],
          breaks = seq(min(biopics.NA$year),
                       max(biopics.NA$year) + (5 - max(biopics.NA$year) %% 5)
                       ,5))

spineMiss(biopics.NA[, c('year','race')],
          breaks = seq(min(biopics.NA$year),
                       max(biopics.NA$year) + (5 - max(biopics.NA$year) %% 5)
                       ,5))

# Country also shows some notable results
spineMiss(biopics.NA[, c('country','box.office')])

# Same goes for color and race
spineMiss(biopics.NA[, c('color','race')])
# This may be related to the year
marginmatrix(biopics.NA[, c('year','race','color','country')], alpha=0.6)

# Use bloxplots to look at the triple relation between year, race and box office
pbox(biopics.NA[,c('year','race','box.office')], selection="none")

# Coordinate plot
parcoordMiss(biopics.NA[,c('race','color','box.office','year')])
             
# Multivariate matrix plot
scattmatrixMiss(biopics.NA[,c('year','race','box.office')], alpha = 0.6)

# Perform single imputation via kNN, to plot imputed values
k = 5
biopics.single.kNN <- kNN(biopics.NA, k = k, imp_var = TRUE)
marginmatrix(biopics.single.kNN[, c(c('year','race','box.office'),
                                    c('year_imp','race_imp','box.office_imp'))], alpha=0.6,
             delimiter = '_imp')
# Intersting: for larger k, box.offices has generally higher imputed values

################################################################################
## Evaluate in-sample performance                                              #
################################################################################

## Run regression with dropped NA observations
# Get regression formula
formula <- paste0('box.office~',
                  paste(names(biopics.NA[,!(names(biopics.NA) %in% c('box.office','title'))]),
                        collapse='+'))
fit.droppedNA <- lm(formula,data = na.omit(biopics.NA))
summary(fit.droppedNA)

## Perform single imputation via kNN and use bootstrap for significance
set.seed(201704241)
k = 5
fit.kNN <- bootstrap(x = biopics.NA[,inc.vars], y.name = 'box.office', 
                     R = 500, k = k, method = 'MM')

## Perform multiple, iterative model-based imputation
set.seed(2017040)
xList <- multimp(x=biopics.NA[,inc.vars], m = 20) 
# use higher m, to compensate for potential lack of convergence
fitList <- fit(xList = xList, y.name = 'box.office', method = 'MM')
fit.MultImp <- pool(fitList, return.more = TRUE, method = 'MM')
# Check how many sets are used
fit.MultImp$imp.sets
fit.MultImp$output.table


# Inspect robustness weights 
rweights <- fit.MultImp$mean.rweights
table(cut(rweights,breaks = seq(0,1,0.1)))
# Being very conservative (choose threshold = 0.5), we get roughly 12 potential outliers 

outlier <- ifelse(rweights<0.5,'red','gray')
miss <- ifelse(is.na(biopics.NA$box.office),'red','gray')
plot(biopics.single.kNN[,c('box.office','year')],col = outlier)
plot(biopics.single.kNN[,c('box.office','year')],col = miss)

################################################################################
## Examine predictive performance                                              #
################################################################################

## Note:
# We are to construct many random training and test sets
# For all methods, imputation is first performed on the training set only
# The test set is imputed with kNN for all methods

# First, inspect how many observations can be regarded as outlier

OutOfSamplePredictions <- function(x, y.name, k = 5, R = 100, method = 'MM', 
                                   training.sample = 0.7, iterations = 100) {
    
    # LTS gives weight 1 to 257 (90.2%) observation in the reweighted estimator. Hence around 10 percent of the observation are outliers.
    # Therefore we compute the 90 percent quantile of the absolute in-sample prediction errors of the MM regression
    # and use this value as tuning parameter c for the tukey bisquared loss function, which is used to determine the out-of-sample prediction errors

    median.error <- as.data.frame(matrix(0, iterations, 3))
    names(median.error) <- c("DroppedNA", "SiBootstrap", "MultImp")
    prediction.errors    <- list(DroppedNA = NULL, SiBootstrap = NULL, MultImp = NULL)
    n.obs <- NROW(x)
    n.train <- ceiling(training.sample * n.obs)
    formula <- paste0('box.office~',
                      paste(names(x[,!(names(x) %in% y.name)]),
                            collapse='+'))
    
    for (i in 1:iterations) {
        #Randomly permute rows of dataset.seed(i)
        x.perm <- x[sample(n.obs),]
        x.train <- x.perm[1:n.train,]

        # Note that it is somewhat easier to directly filter only explanatory
        # Variables from x

        # Model to model.matrix
        x.test  <- x.perm[(n.train+1):n.obs,] 
        x.test <- model.matrix(~., data = x.test)
        y.test <- x.test[,y.name] 
        # Easier to drop response vatiable
        x.test <- x.test[,(colnames(x.test) != y.name)] 

        
        ## Run models
        # With dropped NA
        fit.droppedNA <- lmrob(formula,data = na.omit(x.train), method = 'MM',k.max = 200)
        
        ## Perform single imputation via kNN and use bootstrap for significance
        fit.kNN <- bootstrap(x = x.train, y.name = y.name, R = R, k = k, method = method)
        
        ## Perform multiple, iterative model-based imputation
        xList <- multimp(x = x.train) # Uses default m as percentage missing
        fitList <- fit(xList = xList, y.name = 'box.office', method = 'MM')
        fit.MultImp <- pool(fitList, return.avg.weights.YN = FALSE) 
        
        ## Get prediction error on y.test
        # Problem with type medicine, only one observation with no other missings
        # Solution: perform prediction only on variables present in training set
        names(fit.droppedNA$coefficients)
        prediction.errors$DroppedNA <- y.test - 
            x.test[,names(fit.droppedNA$coefficients)] %*% fit.droppedNA$coefficients
        prediction.errors$SiBootstrap <- y.test - x.test %*% fit.kNN[,1]
        prediction.errors$MultImp <- y.test - x.test %*% fit.MultImp[,1]
        
        median.error[i,] <- sapply(prediction.errors, 
                                   function(x){median(abs(x)) })
        noquote(paste0('iteration: ',i))
    }
    return(median.error)
}

set.seed(2017042413)
pred.error <- OutOfSamplePredictions(x = biopics.NA[,inc.vars], 
                                     y.name = 'box.office', 
                                     k = 5, 
                                     R = 300, 
                                     method = 'MM', 
                                     training.sample = 0.7, 
                                     iterations = 100)

# Get average of median absolute prediction error
colMeans(pred.error)