# Pool point estimates and standard errors
# Input:
# fitList  list of fitted regression models as returned by function fit()
# Output:
# A matrix that contains the pooled coefficients in the first column, the
# pooled standard errors in the second column, the t-statistic in the third
# column, the estimated degrees of freedom in the fourth column, and the
# p-value in the fifth column (see slide 50 of Lecture 5 for an example)
pool <- function(fitList, return.avg.weights.YN = FALSE,...) {
    
    # Retrieve useful parameters
    names.var  <- names(fitList[[1]][['coefficients']])
    n.obs <- NROW(fitList[[1]]$x)
    n.var  <- length(fitList[[1]]$coefficients)
    m <- length(fitList)
    
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
    if (return.avg.weights.YN == TRUE){
        mean.rweights <- sapply(seq(1:n.obs), function(i)
        {
            rweights <- sapply(fitList, function(x){x[['rweights']][i]})
            mean(rweights)
        })        
        output <- list(output.table = output.table,
                       mean.rweights = mean.rweights)
    }else{
        output <- output.table
    }
    
    return(output)
}