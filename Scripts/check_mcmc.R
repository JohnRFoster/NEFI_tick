library(ecoforecastR)
library(runjags)

source("Functions/split_out.R")

type <- "ind"
# type <- "hb"

start <- 1
num.chains <- c(1, 2, 3, 4, 5)

force <- TRUE
cat("force =", force, "\n")

gbr.thresh <- 1.06        # convergence threshold
min.effect.size <- 4000   # minimum effective sample size threshold

if(type == "ind"){
  ## Independent Models ##
  sites <- c("Green", "Henry", "Tea")
  
  # read array job number to paste into output file
  xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
  site <- sites[xx]

  dir <- paste0("../FinalOut/A_Correct/Mice/NymphToAdult/", site)
  model <- paste0("MiceNymphToAdult_", site, "Control")
  # model.save <- paste0("RhoStartMonthEffect_LifeStageDiffPriors_", site, "Control")
  model.save <- model

} else if(type == "hb"){
  
  dir <- "../FinalOut/A_Correct/Mice/LarvaToNymph/"
  model <- "NymphVPD_HB_K_set"
  model.save <- model
  continue <- TRUE

}
cat("Running mcmc diagnostics on", model, "\n")


files <- list.files(dir)

# check if already completed
finished <- grepl("Combined", files)
if(any(finished)) stop("Already combined!", call. = FALSE)


# count number of segments per chain
split.chains <- rep(NA, length(num.chains))
for(c in num.chains){
  pattern <- paste0(model, c)
  split.chains[c] <- length(grep(pattern, files))
}

# check if there are enough chains to combine
num.out <- min(split.chains)
cat("Number of out segments:", num.out, "\n")
if(num.out <= 1) stop("Not enough segments to combine for at least one chain!", call. = FALSE)

load.dir <- file.path(dir, model)

iter.run <- 100000 # number of iterations in each 'out' segment 
total.iter <- iter.run * num.out

cat("Total iterations in each chain:", total.iter, "\n")
cat("Chains being combined:", num.chains, "\n")
cat("Convergence threshold:", gbr.thresh, "\n")
cat("Minimum effective sample size:", min.effect.size, "\n")

spacer <- "_" # for tick models

params <- list()
for(c in seq_along(num.chains)){
  c.load <- num.chains[c]
  load <- paste0(load.dir, c.load, spacer, start, ".RData")
  load(load)
  # if(continue){
  #   out <- split_out(jags.out)
  # }
  chain <- out$out$params
  for(i in (start+1):num.out){
    load <- paste0(load.dir, c.load, spacer, i, ".RData")
    load(load)
    # if(continue){
    #   out <- split_out(jags.out)
    # }
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
if(is.na(burnin) & !force){
  stop("Model not converged!", call. = FALSE)
} else {
  cat("Burnin after:", burnin, "iterations\n")  
}

## determine the segment that burnin occurred
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
  
  file <- paste0("Combined_BURN_", model.save, ".RData")
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
  
  file <- paste0("Combined_thinMat_", model.save, ".RData") 
  save(params.mat, predict.mat, # m.mat,
       file = file.path(dir, file))
  cat("Burned-in thinned matrix saved\n")
  
  cat("Calculating summary statistics\n")
  summary(params.burn)
}

cat("\n--- DONE ---\n")
