print("                                                                 ")
print("-----------------------------------------------------------------")
print("                         Start SCRIPT                            ")
print("-----------------------------------------------------------------")

library(ecoforecastR)

#### model ####
source("Models/GDD_Threshold_HB_K_Window.R") 

## model options
met.variable <- NULL          # met driver on survival
n.adapt <- 150000                  # adaptive iterations
n.chains <- 1                     # number of chains
n.iter <- 100000                  # number of iterations post burnin
iter2save <- 1000                 # number of iterations to save 
thin <- round(n.iter/iter2save)   # thinning interval

## file path to output folder
out.folder <- "../FinalOut/HB_Partial_GDD/K_estimate/WindowLoop"

## model name
out.name <- "GDDSwitch_Window_k_c"

## create file path for output
out.path <- file.path(out.folder,out.name)

## compile model
return <- run_model(n.adapt,n.chains)

## read array job number to paste into output file
chain.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) 

for(samp in 1:500){
  
  jags.out <- coda.samples(model = return$j.model,
                           variable.names = return$monitor,
                           thin = thin,
                           n.iter = n.iter)
  
  time <- Sys.time() - return$start.time 
  cat("\n", "Chain number:", chain.num, "\n", n.iter*samp, "iterations complete after\n")
  print(time)
  
  ## split output
  out <- list(params = NULL, predict = NULL, m.cols = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  m.cols <- grep("m[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$m.cols <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, m.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -c(m.cols, pred.cols)])
  out <- list(out = out, j.model = return$j.model)
  
  iter <- paste(chain.num, samp, sep = "_")
  
  save(out, file = paste(out.path, iter, ".RData", sep = ""))

}

print("-----------------------------------------------------------------")
print("                           END SCRIPT                            ")
print("-----------------------------------------------------------------")
print("                                                                 ")