
convergence_check <- function(jags.out, model, 
                       n.iter = 10000, min.eff.size = 5000, GBR.thresh = 1.02){
  
  enough.samples <- converge <- FALSE
  
  ## split output
  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
  
  # convergence check on parameters
  cat("Calculating PSRF\n")
  GBR.vals <- gelman.diag(out$params, multivariate = FALSE)
  GBR.vals
  converge <- max(GBR.vals$psrf) < GBR.thresh
  cat("Convergence:", converge, "\n")
  
  if(converge){
    cat("Determining burnin\n")
    GBR <- gelman.plot(out$params)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>GBR.thresh,1,any)),1)+1]
    if(is.na(burnin) | burnin > 7500){
      cat("Model not converged!\n")
    } else {
      cat("Burnin after:", burnin, "iterations\n")  
      out$params <- window(out$params, start = burnin)
      out$predict <- window(out$predict, start = burnin)
      enough.samples <- min(effectiveSize(out$params)) >= min.eff.size
      cat("Enough samples:", enough.samples, "\n")
    }
  }
  
  counter <- 1
  while(!converge & !enough.samples){
    counter <- counter + 1
    cat("coda.samples call number:", counter, "\n")
    jags.out <- coda.samples(model = model$j.model,
                             variable.names = model$monitor,
                             n.iter = 10000)
    
    cat("coda samples done, checking mcmc \n")
    
    ## split output
    out <- list(params = NULL, predict = NULL)
    mfit <- as.matrix(jags.out, chains = TRUE)
    pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
    chain.col <- which(colnames(mfit) == "CHAIN")
    out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
    out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
    
    # convergence check on parameters
    cat("Calculating PSRF\n")
    GBR.vals <- gelman.diag(out$params, multivariate = FALSE)
    GBR.vals
    converge <- max(GBR.vals$psrf) < GBR.thresh
    cat("Convergence:", converge, "\n")
    if(!converge) next
    
    cat("Determining burnin\n")
    GBR <- gelman.plot(out$params)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>GBR.thresh,1,any)),1)+1]
    if(is.na(burnin) | burnin > 7500){
      cat("Model not converged!\n")
    } else {
      cat("Burnin after:", burnin, "iterations\n")  
      out$params <- window(out$params, start = burnin)
      out$predict <- window(out$predict, start = burnin)
      enough.samples <- min(effectiveSize(out$params)) >= min.eff.size
      cat("Enough samples:", enough.samples, "\n")
    }
  }
  return(out)
}


