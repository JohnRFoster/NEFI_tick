library(ecoforecastR)
library(runjags)

## Independent Models ##
sites <- c("Green", "Henry", "Tea")
num.out <- c(10, 9, 20)

# read array job number to paste into output file
xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) 

site <- sites[xx]
num.out <- num.out[xx]

dir <- paste0("../FinalOut/A_Correct/ObsModel/Obs.Obs.ObsVert/", site)
model <- paste0("Obs_L1.N1.A2vert_", site, "Control")

# dir <- paste0("../FinalOut/HB_Obs1_Proc1/VPDProc")
# model <- paste0("VPD_ObsProc_beta_111_K_set")

cat("Running mcmc diagnostics on", model, "\n")

load.dir <- file.path(dir, model)

iter.run <- 100000 # number of iterations in each 'out' segment 
num.chains <- c(1, 2, 3, 4, 5)


gbr.thresh <- 1.01
min.effect.size <- 4000

force <- FALSE
cat("force =", force, "\n")

total.iter <- iter.run * num.out
cat("Total iterations in each chain:", total.iter, "\n")
cat("Chains being combined:", num.chains, "\n")
cat("Number of out segments:", num.out, "\n")
cat("Convergence threshold:", gbr.thresh, "\n")
cat("Minimum effective sample size:", min.effect.size, "\n")

spacer <- "_" # for tick models

params <- list()
for(c in seq_along(num.chains)){
  c.load <- num.chains[c]
  load <- paste0(load.dir, c.load, spacer, 1, ".RData")
  load(load)
  chain <- out$out$params
  for(i in 2:num.out){
    load <- paste0(load.dir, c.load, spacer, i, ".RData")
    load(load)
    chain <- combine.mcmc(mcmc.objects = list(chain, out$out$params))
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




