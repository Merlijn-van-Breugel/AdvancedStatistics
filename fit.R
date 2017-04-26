# Fit regression models
# Input:
# xList .. list of imputed data sets as returned by function multimp()
# y.name . name of column to use as y variable 
# method . character string specifying whether to use OLS or MM-estimator
# ........ additional arguments to be passed to modeling function (for example,
# ........ control parameters for the MM-algorithm)
# Output:
# A list of regression models (each fitted to one of the imputed data sets)
fit <- function(xList, y.name = 'V1', method = 'MM', method.inputList = NA, ...){
    
    # Split y and x from xList
    yList     <- lapply(xList, function(x){x[names(x)==y.name]})
    
    ## Perform selected regression method
    # Create formula for regression
    fm <- paste0(y.name,'~',paste0(names(xList[[1]])[(names(xList[[1]]) != y.name)],collapse = '+'))
    # Fit regression models
    if (method == 'MM'){
        fit <- lapply(xList, function(x){lmrob(fm,data = x, method = 'MM',setting = 'KS2014')})
    }else if (method == 'OLS'){
        fit <- lapply(xList, function(x){lm(fm,data = x)})
    }
    
    return(fit) 
}