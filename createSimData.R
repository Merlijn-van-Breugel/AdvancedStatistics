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
                y[outlying.obs]             <- rnorm(prop.outlier*n.obs, 
                                                     mean = 0.1*mean.y, 
                                                     sd = 2)
            }
            outlying.YN[outlying.obs,outlying.var] <- 1
            n.outliers <- sum(outlying.YN)
        }
    }
    
    if (!(prop.miss == 0 | mechanism.miss == 'None')){
        ## Replace prop.miss by NA
        if (mechanism.miss == 'MCAR'){
            # Draw random samples from [0,1] uniform distribution
            mcar <- matrix(runif(size, min=0, max=1),n.obs,n.var)
            # All values smaller than the proportion of missings are replaced by NA
            miss <- ifelse(mcar>prop.miss,1,NA)
            x.obs <- x * miss
        }else if(mechanism.miss == 'MAR'){
            # Create MAR as follows: a covariate in x is missing if y is greater than its miss.prop quantile
            y.quantiles <- quantile(y, seq(from = .1, to = .9, by = .1)) 
            miss <- (y >= y.quantiles[10-floor(n.var*prop.miss*10)])
            x.obs <- x
            for (i in 1:n.obs){
                row.miss <- sample(1:n.var,1)
                x.obs[which(miss)[i],row.miss] <- NA 
            }
        }else if(mechanism.miss == 'MNAR'){
            # Values are missing conditional on their own value
            x.quantiles <- apply(x, 2, function(x){
                quantile(x, seq(from = .1, to = .9, by = .1))})
            x.obs <- x
            for (i in 1:n.obs){
                miss <- ifelse(x.obs[i,] > x.quantiles[10-floor(prop.miss*10),],NA,1)
                x.obs[i,] <- x.obs[i,] * miss
            }
        }
    }else{x.obs <- x}
    sum(is.na(x.obs)) #This should be roughly equal to n.obs * n.var * prop.miss
    return(list(
        data = cbind(y,x.obs),
        prop.miss = sum(is.na(x.obs))/size))
}