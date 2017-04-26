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
            fit <- lmrob(fm,data = sample.imp, method = 'MM')
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
    t.stat <- coef.mean/sd
    df = n.obs - n.var - 1
    # Two-sided, get p value
    p.val <- 2*pt(abs(t.stat), df = df, lower = FALSE) 
    # Nice function to get significance
    symp <- symnum(p.val, corr = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))
    output.table <- (data.frame('Estimate' = coef.mean,
                                'Std.Error' = sd,
                                't.value' = t.stat,
                                'p.val' = p.val,
                                'sign' = symp))
    return(output.table)    
}