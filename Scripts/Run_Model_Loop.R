print("                                                                 ")
print("-----------------------------------------------------------------")
print("                         Start SCRIPT                            ")
print("-----------------------------------------------------------------")

library(ecoforecastR)

#### model ####
source("Models/GDD_K_set.R") 

## model options
site.run <- "Tea Control"
met.proc <- NULL                 # met driver on survival
n.adapt <- 300000                 # adaptive iterations
n.chains <- 1                     # number of chains
n.iter <- 100000                  # number of iterations post burnin
iter2save <- 5000                 # number of iterations to save 
thin <- round(n.iter/iter2save)   # thinning interval


## file path to output folder
# out.folder <- "../FinalOut/HB_Obs1_Proc1/ReproMax0"
out.folder <- "../FinalOut/A_Correct/NULL/Tea"

## model name HB Runs
# out.name <- "MaxTemp_ObsProc_beta_111_ReproMax0_K_set"

## Model name Individual site runs
out.name <- paste("NULL", gsub(" ","",site.run),sep="_")

## create file path for output
out.path <- file.path(out.folder,out.name)

## compile model
# return <- run_model(met.proc,n.adapt,n.chains)
# return <- run_model(site.run, n.adapt,n.chains) # for independent fits
return <- run_model(site.run, met.proc, n.adapt,n.chains) # for independent fits with met.proc

## read array job number to paste into output file
chain.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) 

for(samp in 1:5000){
  # for(samp in 1:5){ # for testing
  loop.start <- Sys.time()
  
  jags.out <- coda.samples(model = return$j.model,
                           variable.names = return$monitor,
                           # thin = thin,
                           n.iter = n.iter)
  
  loop.time <- Sys.time() - loop.start
  cat(n.iter, "iterations in\n")
  print(loop.time)
  
  time <- Sys.time() - return$start.time 
  cat("\nChain number:", chain.num, "\n", n.iter*samp, "iterations complete after\n")
  print(time)
  
  ## split output
  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
  out <- list(out = out, j.model = return$j.model)
  
  iter <- paste(chain.num, samp, sep = "_")
  
  variable.names <- return$monitor
  
  save(out, 
       variable.names,
       file = paste(out.path, iter, ".RData", sep = ""))
  
}

print("-----------------------------------------------------------------")
print("                           END SCRIPT                            ")
print("-----------------------------------------------------------------")
print("                                                                 ")