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