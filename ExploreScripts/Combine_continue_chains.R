library(ecoforecastR)

source("Functions/combine_chains.R")
source("Functions/plot_posterior_chains.R")



path <- "../FinalOut/Independent_Fits/GDDThreshold/Window/GDDSwitch_GreenControl"
out <- combine_chains(path)

path <- "../FinalOut/Independent_Fits/GDDThreshold/Window/GDDSwitch_GreenControl_Cont1_"


out.cont <- combine_chains(path)
plot(out.cont$params)

out.win <- window(out$params,start=1,end=nrow(out$params[[1]]),length.out=500)
str(out.win)

seq <- round(seq(1,10000,length.out = 500))
out.win <- list()
out.win <- out$params[[1]][seq,]
out.win <- out$params[[2]][seq,]
out.win <- out$params[[3]][seq,]
out.win <- out$params[[4]][seq,]
out.win <- out$params[[5]][seq,]

out.all <- as.mcmc(as.mcmc(out),as.mcmc(out.cont))

out1 <- out$params
out2 <- out.cont$params

c <- list()
c[[1]] <- as.mcmc(rbind(out$params[[1]],out.cont$params[[1]]))
c[[2]] <- as.mcmc(rbind(out$params[[2]],out.cont$params[[2]]))
c[[3]] <- as.mcmc(rbind(out$params[[3]],out.cont$params[[3]]))
c[[4]] <- as.mcmc(rbind(out$params[[4]],out.cont$params[[4]]))
c[[5]] <- as.mcmc(rbind(out$params[[5]],out.cont$params[[5]]))

out.all <- as.mcmc(c)
plot(out.all)

effectiveSize(out.all)

ammend_chains <- function(dir, out1, out2, save = NULL){
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

