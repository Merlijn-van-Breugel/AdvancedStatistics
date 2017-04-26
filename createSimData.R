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
                                                     mean = 0.1*mean.y, 
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
        }else if(mechanisms.miss == 'MAR'){
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