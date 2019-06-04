library(ecoforecastR)
library(runjags)

source("Functions/combine_chains.R")
source("Functions/plot_posterior_chains.R")

path <- "FinalOut/Independent_Fits/GDDThreshold/Window/Temp_GDDSwitch_GreenControl_Cont1_"
chain.numbs <- 1:5

out.cont <- combine_chains(path,chain.numbs)
plot(out.cont$params)


path <- "FinalOut/Independent_Fits/GDDThreshold/Window/Temp_GDDSwitch_GreenControl"
out <- combine_chains(path,chain.numbs)

out.all <- as.mcmc(as.mcmc(out),as.mcmc(out.cont))
out1 <- out$params
out2 <- out.cont$params
c <- list()
c[[1]] <- as.mcmc.list(out$params[[1]],out.cont$params[[1]])
c[[2]] <- as.mcmc.list(out$params[[2]],out.cont$params[[2]])
c[[3]] <- as.mcmc.list(out$params[[3]],out.cont$params[[3]])
c[[4]] <- as.mcmc.list(out$params[[4]],out.cont$params[[4]])
c[[5]] <- as.mcmc.list(out$params[[5]],out.cont$params[[5]])

out.all <- as.mcmc.list(c1,c2,c3,c4,c5)
plot(out.all)
