createBaseData <- function(n.obs = 100, n.var = 3, 
                           mean.x = c(5,10,15), 
                           sd.x = c(1,2,3), sd.y = 2,
                           c = -5, coef = c(1,2,3)){
    
    ## Create base data set
    # Create p-1 explanatory variables that are normally distributed 
    x <- mapply(function(x,y){rnorm(x,y,n=n.obs)},x=mean.x,y=sd.x)
    colnames(x) <- paste0('x',1:n.var)
    # Compute y, given the pre-defined coefficients
    y <- rnorm(n, mean = c + x %*% coef, sd = sd.y)
    
    return(cbind(y,x))
}