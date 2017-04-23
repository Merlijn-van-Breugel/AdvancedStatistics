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
install.packages('doParallel')  
library(doParallel)  
require(snow)
require(VIM)
require(robustbase)   
# Load data
load("biopics.RData")
# Set seed. Best-practice: large seed for more randomness
set.seed(42424242)


################################################################################
## Create functions                                                            #                 
################################################################################

## Function to create dataset with missing data and outliers
createData <- function()
# Create data set with missings
n    <- 100 
ncol <- 5
x <- matrix(rexp(n*ncol, rate=.1), ncol=ncol)
# Drop values
miss.perc <- 10
n.miss <- NROW(x)
x <- as.data.frame(apply(x, 2, function(x) {x[sample(c(1:n.miss), floor(n.miss/ncol))] <- NA; x} ))
x


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
multimp <- function(x, m = NULL, mixed = NULL, robust = FALSE, robMethod = 'MM', init.method = 'kNN',
                    ...) {
    
    # Set a sensible default for the number of imputations m
    # Follow advice of Bodner (2008) and White et al. (2011), set m default to.. 
    # ... percentage of missing information
    if (is.null(m)){
        m <- 5*round(100*mean(is.na(x))/5)
    }
    ## Use function irmi() from VIM package multiple imputations
    
    # Check whether there are semi-continuous variables
    # Here we use a rather ambiguous threshold, but this is more of a sanity check
    if (is.null(mixed)){
        mixed <- apply(x, 
                    2, 
                    function(x) {(sum(x == 0, na.rm = TRUE)/NROW(x)) > 0.2 } )
        if (sum(mixed) == 0){
            mixed <- NULL}
    }   
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

xList <- multimp(x)

# Fit regression models
# Input:
# xList .. list of imputed data sets as returned by function multimp()
# y.name . name of column to use as y variable 
# method . character string specifying whether to use OLS or MM-estimator
# ........ additional arguments to be passed to modeling function (for example,
# ........ control parameters for the MM-algorithm)
# Output:
# A list of regression models (each fitted to one of the imputed data sets)
fit <- function(xList, y.name = 'V1', method = 'MM', method.inputList = NA, ...) {
    
    # Split y and x from xList
    yList     <- lapply(xList, function(x){x[names(x)==y.name]})
    
    ## Perform selected regression method
    # Check whether input method is allowed
    if (!(method %in% c('OLS','MM'))){
        warning('Chosen regression method is not supported, MM is chosen instead')
        method <- 'MM'
    }
    # Create formula for regression
    fm <- paste0(y.name,'~',paste0(names(xList[[1]])[(names(xList[[1]]) != y.name)],collapse = '+'))
    # Fit regression models
    if (method == 'MM'){
        fit <- lapply(xList, function(x){lmrob(fm,data = x, method = 'MM')})
    }else if (method == 'OLS'){
        fit <- lapply(xList, function(x){lm(fm,data = x)})
    }

    return(fit) 
}

fitList <- fit(xList = xList, y.name = 'V1')

# Pool point estimates and standard errors
# Input:
# fitList  list of fitted regression models as returned by function fit()
# Output:
# A matrix that contains the pooled coefficients in the first column, the
# pooled standard errors in the second column, the t-statistic in the third
# column, the estimated degrees of freedom in the fourth column, and the
# p-value in the fifth column (see slide 50 of Lecture 5 for an example)
pool <- function(fitList, ...) {
    
    # Check whether input only contains OLS and MM models
    log <- sapply(fitList, function(x){
                    method <- (summary(x))[['call']][['method']]
                    if (!(method %in% c('OLS','MM'))){
                    stop('Regression model input has to be OLS or MM model')
                    }
    })
    # Retrieve useful parameters
    names.var  <- names(fitList[[1]][['coefficients']])
    n.obs <- NROW(fitList[[1]]$x)
    n.var  <- length(fitList[[1]]$coefficients)
    m <- length(fitList)
    
    # Start with pooled coefficients estimates
    mean.pool <- sapply(seq(1:n.var), function(i)
    {
        # m <- sapply(lapply(fitList, `[[`, 'coefficients'), function(x) x[i])
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
                                - coef.mean[i])^2})
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
    t.stat <- coef.mean/sd
    # Two-sided, get p value
    p.val <- 2*pt(abs(t.stat), df = df, lower = FALSE) 
    # Nice function to get significance
    symp <- symnum(p.val, corr = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))
    output.table <- as.matrix(data.frame('Estimate' = coef.mean,
                               'Std.Error' = sd,
                               't.value' = t.stat,
                               'df' = df,
                               'p.val' = p.val))
    rownames(output.table) <- names.var
    return(output.table)
}

output <- pool(fitList)

## Function for the bootstrap with kNN imputation and linear regression
# Input:
# x       data set with missing values
# R       number of replications
# k       number of neighbors for kNN imputation
# method  character string specifying whether to use OLS or MM-estimator
# ...     additional arguments to be passed to modeling function (for example,
#         control parameters for the MM-algorithm)
# Output:
# A list with the following components:
# replicates  a matrix with all coefficient estimates from all replications
# summary     a matrix that contains the point estimates of the coefficients in
#             the first column, the standard errors in the second column, the
#             z-statistic in the third column, and the p-value in the fourth
#             column (see slide 29 of Lecture 5 for an example)
bootstrap <- function(x, y.name = 'V1', R = 1000, k = 5, method = 'MM', ...) {
    
    # Check whether input method is allowed
    if (!(method %in% c('OLS','MM'))){
        warning('Chosen regression method is not supported, MM is chosen instead')
        method <- 'MM'
    }
    
    # Retrieve useful parameters
    names.var  <- names(x)
    n.obs <- NROW(x)
    n.var <- (NCOL(x)-1)
    # Already compute formula for regression
    fm <- paste0(y.name,'~',paste0(names(x)[(names(x) != y.name)],collapse = '+'))
    
    ## Get R bootstrap samples with regression results
    
    # Run in parallel for drastic performance increase
    cores <- detectCores()
    registerDoParallel(detectCores())  
    # Clusters willed be stopped automatically
    
    # time.start <- Sys.time()
    bootstrap <-  foreach(i=1:R, .combine=rbind) %dopar% {  
        # Load packages into parallel workers. Maybe there is a more elegant way to do this
        if (i %in% 1:cores){
            library(VIM)
            library(robustbase)
        }
        
        # Sample n observations with replacement
        s <- x[sample(NROW(x), n.obs, replace = TRUE),]
        # Impute sample to get complete data matrix
        
        # We use kNN method from the VIM package for single imputation
        sample.imp <- kNN(s, k = k, imp_var = FALSE)
        
        # Fit regression models
        if (method == 'MM'){
            fit <- lmrob(fm,data = sample.imp, method = 'MM')
        }else if (method == 'OLS'){
            fit <- lm(fm,data = sample.imp)
        }
        fit[['coefficients']] 
    }
    # time.end <- Sys.time()
    # time.end - time.start 
    
    coef.mean <- colMeans(bootstrap)

    # Compute standard deviation of bootstrap replicates
    sd <- sqrt((1/(R-1))*colSums(sweep(bootstrap,2,coef.mean)^2))  

    # Get t-statistics based on pooled mean and variance
    t.stat <- coef.mean/sd
    df = n.obs - n.var - 1
    # Two-sided, get p value
    p.val <- 2*pt(abs(t.stat), df = df, lower = FALSE) 
    output.table <- as.matrix(data.frame('Estimate' = coef.mean,
                                         'Std.Error' = sd,
                                         't.value' = t.stat,
                                         'p.val' = p.val))
    rownames(output.table) <- names.var
    return(output.table)    
}

test <- bootstrap(x, k = 5, R = 10)
