#' This takes the continued jags run and ammends it to the previous run
#'
#'@param dir file path to model output
#'@param chain.numbs vector of chain numbers

ammend_chains <- function(dir, out1, out2, thin = 10000, save = NULL){
  # storage
  params <- state <- state.m <- list()
  
  # JAGS output
  part1 <- combine_chains(file.path(dir,out1))
  part2 <- combine_chains(file.path(dir,out2))
  
  # number of chains
  n.chains <- length(part1$params)
  
  for(c in 1:nchains){
    params[[c]] <- as.mcmc(rbind(part1$params[[c]],part2$params[[c]]))
    state[[c]] <- as.mcmc(rbind(part1$predict[[c]],part2$predict[[c]]))
    state.m[[c]] <- as.mcmc(rbind(part1$predict.m[[c]],part2$predict.m[[c]]))
  }
  
  out <- list(params = as.mcmc(params),
              predict = as.mcmc(state),
              predict.m = as.mcmc(state.m))
  
  if(!is.null(save)){
    save(out,file.path(dir,save))
  }
  return(out)
}