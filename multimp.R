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