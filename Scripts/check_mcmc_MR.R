library(ecoforecastR)
library(runjags)


site <- "Tea"

# dir <- paste0("../FinalOut/Independent_Fits/GDDThreshold/RH_ObsProc/beta_111/", site)
# model <- paste0("MaxRH_ObsProc_beta_111_K_set_", site, "Control")

dir <- paste0("../FinalOut/HenryControlMR")
model <- "HenryControlMR_Recruit_"

cat("Running mcmc diagnostics on", model, "\n")

load.dir <- file.path(dir, model)

iter.run <- 1000 # number of iterations in each 'out' segment 
num.chains <- c(1, 2, 3, 4, 5)

# 1 = 80
# 2 = 56
# 3 = 59
# 4 = 42
# 5 = 59

num.out <- 3
gbr.thresh <- 1.04
min.effect.size <- 4000

force <- FALSE
cat("force =", force, "\n")

total.iter <- iter.run * num.out
cat("Total iterations in each chain:", total.iter, "\n")
cat("Chains being combined:", num.chains, "\n")
cat("Number of out segments:", num.out, "\n")
cat("Convergence threshold:", gbr.thresh, "\n")
cat("Minimum effective sample size:", min.effect.size, "\n")

spacer <- "." # for mice models

## split output function
split_out <- function(jags.out){
  out <- list(params = NULL, other = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  gamma.cols <- grep("gamma[", colnames(mfit), fixed = TRUE)
  lambda.cols <- grep("lambda", colnames(mfit), fixed = TRUE)
  theta.cols <- grep("theta", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col,gamma.cols,lambda.cols,theta.cols)])
  out$other <- ecoforecastR::mat2mcmc.list(mfit[, -c(gamma.cols,lambda.cols,theta.cols)])
  out <- list(out = out)
  return(out)
}

params <- list()
for(c in seq_along(num.chains)){
  c.load <- num.chains[c]
  load <- paste0(load.dir, c.load, spacer, 1, ".RData")
  load(load)
  out <- split_out(jags.out)
  chain <- out$params
  for(i in 2:num.out){
    load <- paste0(load.dir, c.load, spacer, i, ".RData")
    load(load)
    chain <- combine.mcmc(mcmc.objects = list(chain, out$params))
    if(i %% 5 == 0){
      cat(i, "segments combined\n")
    }
  }
  params[[c]] <- as.mcmc(chain)
  cat("Chain", c.load, "complete\n")
}

cat("All chains combined\n")
params <- as.mcmc.list(params)

cat("Calculating PSRF\n")
GBR.vals <- gelman.diag(params, multivariate = FALSE)
converge <- max(GBR.vals$psrf) < gbr.thresh
GBR.vals

cat("Determining burnin\n")
GBR <- gelman.plot(params)
burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>gbr.thresh,1,any)),1)+1]
if(is.na(burnin)){
  stop("Model not converged!", call. = FALSE)
} else {
  cat("Burnin after:", burnin, "iterations\n")  
}

## determine the segment that burnin occured
start <- ceiling(burnin / iter.run)
params.burn <- window(params, start = (start-1)*iter.run+1)

cat("Calculating effective sample size\n")
effsize <- effectiveSize(params.burn)
enough.samples <- min(effsize) > min.effect.size
if(enough.samples | force){
  print(effsize)
  if(!enough.samples){
    cat("Not enough samples for at least one parameter\n")
  } else {
    cat("Minimum effective sample sized reached\n")
  }
} else {
  print(effsize)
  cat("\n")
  stop("Not enough samples for at least one parameter", call. = FALSE)
}

## if there is burnin and convegence and enough samples
if((!is.na(burnin) & converge & enough.samples) | force){
  
  # plot trace plots
  # trace <- paste0(model, "_params.jpeg")
  # jpeg(filename = file.path(dir, trace),
  #      width = 1200,
  #      height = 1200,
  #      pointsize = 40,
  #      quality = 100)
  # plot(params.burn)
  # dev.off()
  
  ## combine predict and m.cols
  cat("Combining predict output\n")
  predict <- list()
  m.cols <- list()
  for(c in seq_along(num.chains)){
    c.load <- num.chains[c]
    load <- paste(load.dir, c.load, spacer, start, ".RData", sep = "")
    load(load)
    chain.p <- out$out$predict
    # chain.m <- out$out$m.cols
    for(i in (start + 1):num.out){
      load <- paste(load.dir, c.load, spacer, i, ".RData", sep = "")
      load(load)
      chain.p <- combine.mcmc(mcmc.objects = list(chain.p, out$out$predict))
      # chain.m <- combine.mcmc(mcmc.objects = list(chain.m, out$out$m.cols))
      if(i %% 5 == 0){
        cat("\t", i, "segments combined\n")
      }
    }
    predict[[c]] <- as.mcmc(chain.p)
    # m.cols[[c]] <- as.mcmc(chain.m)
    cat("\t", "Chain", c.load, "complete\n")
  }
  
  ## remove burnin and save output
  
  predict.burn <- as.mcmc.list(predict)
  # m.burn <- as.mcmc.list(m.cols)
  
  file <- paste0("Combined_BURN_", model, ".RData")
  save(params.burn, predict.burn, # m.burn,
       file = file.path(dir, file))
  cat("Burned-in mcmc file saved\n")
  
  params.mat <- as.matrix(params.burn)
  predict.mat <- as.matrix(predict.burn)
  # m.mat <- as.matrix(m.burn)
  
  iter <- nrow(params.mat)
  thin <- seq(1, iter, by = ceiling(iter / 10000))
  
  params.mat <- params.mat[thin,]
  predict.mat <- predict.mat[thin,]
  # m.mat <- m.mat[thin,]
  
  file <- paste0("Combined_thinMat_", model, ".RData") 
  save(params.mat, predict.mat, # m.mat,
       file = file.path(dir, file))
  cat("Burned-in thinned matrix saved\n")
  
  cat("Calculating summary statistics\n")
  summary(params.burn)
}



cat("\n--- DONE ---\n")




