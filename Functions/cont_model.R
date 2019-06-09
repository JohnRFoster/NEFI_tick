

library(ecoforecastR)

cont_model <- function(path, model, site, num.cont){

  chain <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file
  
  site <- gsub(" ","",site)
  
  load <- paste(file.path(path,model),site,chain,".RData",sep="")
  
  load(load)
  
  j.model <- out$j.model
  
  j.model$recompile()
  
  
  burnin <- 0                   # burnin iterations
  n.iter <- 800000                  # number of iterations post burnin
  iter2save <- 10000                 # number of iterations to save 
  thin <- round(n.iter/iter2save)   # thinning interval
  
  monitor <- c("x","m",
               "deviance",
               "phi.l.mu",
               "phi.n.mu",
               "phi.a.mu",
               "grow.ln.mu",
               "grow.na.mu",
               "repro.mu",
               "SIGMA",
               "beta.11",
               "beta.22",
               "beta.33",
               "theta.larvae",
               "theta.nymph",
               "theta.adult")
  
  load.module("glm")
  load.module("dic")
  
  jags.out <- coda.samples(model = j.model,
                           variable.names = monitor,
                           burnin = burnin,
                           thin = thin,
                           n.iter = n.iter)
  
  out <- list(params = NULL, predict = NULL, m.cols = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  m.cols <- grep("m[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$m.cols <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, m.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -c(m.cols, pred.cols)])
  out <- list(out = out, j.model = j.model)

  out.path <- paste(file.path(path,model),site,"_Cont",num.cont,"_",sep="")

  save(out, file = paste(out.path, chain, ".RData", sep = ""))

}

path <- "../FinalOut/Independent_Fits/GDDThreshold/LowThreshOnly"
model <- "Temp_GDDSwitch_gammaBeta33_"
site <- "Henry Control"

cont_model(path, model, site, 1)
