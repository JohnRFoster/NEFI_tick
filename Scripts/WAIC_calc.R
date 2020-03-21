library(ecoforecastR)
library(LaplacesDemon)

source("Functions/site_data_met.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/ua_parts.R")
source("Functions/obs_prob.R")

data.hb <- cary_ticks_met_JAGS()

site.folders <- c("Green", "Henry", "Tea")
sites <- c("Green Control", "Henry Control", "Tea Control")
top.dir <- "../FinalOut/A_Correct/ObsModel/Obs.Obs.ObsVert/"
model.name <- "Combined_thinMat_Obs_L1.N1.A2vert_"
met.variable <- NULL

for(i in 1:3){
  dir <- paste0(top.dir, site.folders[i])
  model <- paste0(model.name, site.folders[i], "Control.RData")
  load(file.path(dir, model))
  
  cat("Calculating WAIC for", model, "\n")

  data <- site_data_met(site = sites[i], 
                        met.variable = met.variable,
                        data = data.hb)
  obs.temp <- data$met.obs
  obs.temp[data$met.obs.miss] <- mean(obs.temp, na.rm = TRUE)
  N_est <- data$N_est
  y <- data$y
  ua <- ua_parts(params.mat, predict.mat, c("parameter", "ic"))
  
  # storage
  like.l <- like.n <- like.a <- matrix(NA, nrow(predict.mat), N_est)
  
  # observation probability
  obs.prob <- obs_prob(ua, N_est, obs.temp)
  num <- length(ua$deviance)
  
  for(t in 1:N_est){ # loop over observations
    
    # zero inflate latent state
    l <- ua$IC[, paste("x[1,",t, "]",sep="")]# * rbinom(num, 1, obs.prob$theta.larva[,t])
    n <- ua$IC[, paste("x[2,",t, "]",sep="")]# * rbinom(num, 1, obs.prob$theta.nymph[,t])
    a <- ua$IC[, paste("x[3,",t, "]",sep="")]# * rbinom(num, 1, obs.prob$theta.adult[,t])
    
    # set zero predictions to 1E-10 (otherwise get NaN in likelihood)
    l[l==0] <- 1E-10
    n[n==0] <- 1E-10
    a[a==0] <- 1E-10
    
    # calculate likelihood
    like.l[,t] <- dgpois(y[1,t], l, 1-obs.prob$theta.larva[,t])
    like.n[,t] <- dgpois(y[2,t], n, 1-obs.prob$theta.nymph[,t])
    like.a[,t] <- dgpois(y[3,t], a, 1-obs.prob$theta.adult[,t])
  }
  
  like <- cbind(like.l, like.n, like.n)
  fbar <- colMeans(like)
  like[like==0] <- 1E-10
  Pw <- sum(apply(log(like),2,var))
  lppd <- sum(log(fbar))
  WAIC <- -2*lppd + 2*Pw
  
  cat("Pw:", Pw, "\n")
  cat("lppd:", lppd, "\n")
  cat("WAIC:", WAIC, "\n\n\n")
  
  waic.dat <- data.frame(site = site.folders[i],
                    Pw = Pw,
                    lppd = lppd,
                    WAIC = WAIC,
                    model = model.name)
  
  load("WAIC.RData")
  dat <- rbind(dat, waic.dat)
  save(dat, file = "WAIC.RData")
  
}

dat.green <- subset(dat, site=="Green")
dat.green$WAIC.diff <- min(dat.green$WAIC) - dat.green$WAIC

dat.henry <- subset(dat, site=="Henry")
dat.henry$WAIC.diff <- min(dat.henry$WAIC) - dat.henry$WAIC

dat.tea <- subset(dat, site=="Tea")
dat.tea$WAIC.diff <- min(dat.tea$WAIC) - dat.tea$WAIC


# dat <- dat[1:36,]
# save(dat, file = "WAIC.RData")








