#' This function combines parallel chains from JAGS 
#'
#'@param path file path to model out
#'@param chain.numbs vector of chain numbers


combine_chains <- function(path, chain.numbs){
  
  # initialize
  params <- predict <- predict.m <- list()
  for(i in 1:length(chain.numbs)){
    # paste model file path to chain number
    chain <- paste(path,chain.numbs[i],".RData",sep="")
    load(chain)
    params[[i]] <- as.mcmc(out$out$params)
    predict[[i]] <- as.mcmc(out$out$predict)
    predict.m[[i]] <- as.mcmc(out$out$m.cols) 
  }
  out <- list(params = as.mcmc(params),
              predict = as.mcmc(predict),
              predict.m = as.mcmc(predict.m))
  return(out)
}


# Used to test - might come in handy later...



# params <- predict <- predict.m <- list()
# load("FinalOut/Independent_Fits/GDDThreshold/Window/Temp_GDDSwitch_HenryControl1.RData")
# params[[1]] <- as.mcmc(out$params)
# predict[[1]] <- as.mcmc(out$predict)
# predict.m[[1]] <- as.mcmc(out$m.cols)
# load("FinalOut/Independent_Fits/WithMolting/Temp_Molting_GreenControl2.RData")
# params[[2]] <- as.mcmc(out$params)
# predict[[2]] <- as.mcmc(out$predict)
# predict.m[[2]] <- as.mcmc(out$m.cols)
# load("FinalOut/Independent_Fits/WithMolting/Temp_Molting_GreenControl3.RData")
# params[[3]] <- as.mcmc(out$params)
# predict[[3]] <- as.mcmc(out$predict)
# predict.m[[3]] <- as.mcmc(out$m.cols)
# load("FinalOut/Independent_Fits/WithMolting/Temp_Molting_GreenControl4.RData")
# params[[4]] <- as.mcmc(out$params)
# predict[[4]] <- as.mcmc(out$predict)
# predict.m[[4]] <- as.mcmc(out$m.cols)
# load("FinalOut/Independent_Fits/WithMolting/Temp_Molting_GreenControl5.RData")
# params[[5]] <- as.mcmc(out$params)
# predict[[5]] <- as.mcmc(out$predict)
# predict.m[[5]] <- as.mcmc(out$m.cols)
# 
# params <- as.mcmc(params)
# predict <- as.mcmc(predict)
# predict.m <- as.mcmc(predict.m)