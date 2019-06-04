library(ecoforecastR)


source("Functions/combine_chains.R")
path <- "FinalOut/Independent_Fits/GDDThreshold"
chain.numbs <- 1:5
#out <- combine_chains(file.path(path,"Temp_GDDSwitch_l2n_lowONLY_GreenControl"),chain.numbs)
out <- combine_chains(file.path(path,"Window/Temp_GDDSwitch_TeaControl"),c(1,2,3,4,5))
plot(out$params)
params <- window(out$params, start = 5000)
plot(params[-5])
effectiveSize(params)
effectiveSize(out$params)

params <- out$params

j.model <- out$j.model
