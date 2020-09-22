library(runjags)
convergence_check <- function(jags.out, model, monitor, n.iter, print = FALSE,
                       min.eff.size = 3000, GBR.thresh = 1.1){
  
  enough.samples <- converge <- FALSE
  
  ## split output
  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
  
  params.summary <- summary(out$params)
  iterations <- params.summary$end
  
  # convergence check on parameters
  cat("Calculating PSRF\n")
  GBR.vals <- gelman.diag(out$params, multivariate = FALSE)
  if(print) print(GBR.vals)
  converge <- max(GBR.vals$psrf) < GBR.thresh
  cat("Convergence:", converge, "\n")
  
  if(converge){
    cat("Determining burnin\n")
    GBR <- gelman.plot(out$params)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>GBR.thresh,1,any)),1)+1]
    if(is.na(burnin) | burnin > iterations*0.75){
      cat("Model not converged!\n")
    } else {
      cat("Burnin after:", burnin, "iterations\n")  
      out$params <- window(out$params, start = burnin)
      out$predict <- window(out$predict, start = burnin)
      enough.samples <- min(effectiveSize(out$params)) >= min.eff.size
      cat("Enough samples:", enough.samples, "\n")
    }
  }
  
  counter <- thin <- 1
  while(!converge | !enough.samples){
    counter <- counter + 1
    cat("coda.samples call number:", counter, "\n")
    
    if(counter %% 10 == 0) { # print every 10 iterations through
      print <- TRUE
      thin <- counter*10 # thinning always increases
    } else {
      print <- FALSE
    }
    
    new.out <- coda.samples(model = model,
                            variable.names = monitor,
                            n.iter = n.iter,
                            thin = thin)
    jags.out <- combine.mcmc(mcmc.objects = list(jags.out, new.out),
                             collapse.chains = FALSE)
    
     
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
    if(print) print(GBR.vals)
    converge <- max(GBR.vals$psrf) < GBR.thresh
    cat("Convergence:", converge, "\n")
    if(!converge) next
    
    cat("Determining burnin\n")
    GBR <- gelman.plot(out$params)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>GBR.thresh,1,any)),1)+1]
    if(is.na(burnin)){
      cat("Model not converged!\n")
    } else {
      cat("Burnin after:", burnin, "iterations\n")  
      out$params <- window(out$params, start = burnin)
      out$predict <- window(out$predict, start = burnin)
      effect.size <- effectiveSize(out$params)
      enough.samples <- min(effect.size) >= min.eff.size
      cat("Enough samples:", enough.samples, "\n")
      # if(!enough.samples & print) print(min(effect.size))
      min.index <- which(effect.size == min(effect.size))
      print(effect.size[min.index])
    }
  }
  return(out)
}


