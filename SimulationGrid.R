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
        length(prop.outlier.range)*
        length(mechanism.miss.range)*
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
    clusterCall(cl, function(){source('createBaseData.R')})
    clusterCall(cl, function(){source('createSimData.R') })
    
    
    # Now loop over all configurations
    for (i in 1:length(mechanism.miss.range)){
        for (j in 1:length(prop.miss.range)){
            for (k in 1:length(mechanism.outlier.range)){
                for (l in 1:length(prop.outlier.range)){
                    # For each configuration, we run this sim times
                    config.run <- foreach(s=1:sim, .combine=rbind) %dopar% {
                        # For the first workers, we have to source our functions and packages
                        # if (s <= detectCores()){
                        #    library(VIM)
                        #    library(robustbase)
                        #    source('createBaseData.R')
                        #    source('createSimData.R')} 
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
                        colnames(sim.data) <- colnames(base.data) # Weird, should not be needed
                        regress.coef <- NULL
                        
                        # Choose k in accordance to convention: sqrt(N)
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
                                                         setting = 'KS2014')
                                    }
                                    coef <- fit.KNN$coefficients
                                }else if (method.imp == 'multiple'){
                                    fitList <- fit(xList = xList, y.name = 'y',method = method.sim)
                                    coef <- pool(fitList)[,1]
                                    if (s == 1){
                                        coef.names <- names(coef)
                                    }
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
    dimnames(res)[[2]] <- colnames(regress.coef)
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