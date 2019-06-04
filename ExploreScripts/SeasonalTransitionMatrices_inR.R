## attempting to use a questing matrix that is aggregated to the daily scale
## and seasonal transition matrices that are used at each observation
## trying to fit this model in R

## this is just the hindcast, using the parameter values from the JAGS historical fits



source("Functions/Future_Met.R")
source("Functions/get_last_day.R")

load("/projectnb/dietzelab/fosterj/FinalOut/HB_Partial_Temp/Out_Temp_Partial_HB.RData")

params <- window(params, start = 3000)
predict <- window(predict, start = 3000)

params.mat <- as.matrix(params)
ic.mat <- as.matrix(predict)

load("/projectnb/dietzelab/fosterj/FinalOut/SurvivalModels/Fed_Larvae_MaxTemp_summer.RData")
molt.nymph <- as.matrix(jags.out)
molt.nymph <- window(molt.nymph, start = 3000)
molt.nymph <- as.matrix(molt.nymph)

forecast_with_molting_temp <- function(type, site, params, ic, N_est, df, Nmc, draw,
                                       q.adult_d.adult,molt.nymph){
  
  # site index
  if(site == "Green Control"){
    s <- 1
  } else if(site == "Henry Control"){
    s <- 2
  } else {
    s <- 3
  }
  
  df <- rep(df, N_est)
  dt.index <- cumsum(df)
  N_days <- dt.index[length(dt.index)]
  
  last <- get_last_day()
  
  n.obs <- last[s,"n.obs"]
  
  met.hind <- future_met(site, 10)
  cum.gdd <- met.hind$cum.gdd[1:N_days]
  temp <- met.hind$temp.scale[1:N_days]
  
  # storage
  pred <- array(dim = c(7,N_est,Nmc))
  A <- array(0, dim = c(3,3,N_days))
  P <- array(0, dim = c(2,2,N_days))
  Q <- array(0, dim = c(2,2,N_days))
  R <- array(0, dim = c(2,2,N_days))
  S <- array(0, dim = c(2,2,N_days))
  
  # mean parameters
  param.mean <- apply(params, 2, mean)
  molt.nymph.surv.mean <- apply(molt.nymph, 2, mean)
  
  # select appropriate initial conditions
  ic <- ic[,grep(paste(n.obs, ",", s, "]", sep = ""), colnames(ic))]
  
  if("ic" %in% type){
    IC <- ic[draw,]
  } else {
    IC.mean <- apply(ic, 2, mean)
    IC.names <- names(IC.mean)
    IC <- matrix(IC.mean, Nmc, ncol = length(IC.mean), byrow = TRUE)
    colnames(IC) <- IC.names
  }
  
  # use appropriate site random effects
  alpha.11 <- params[,paste("alpha.11[", s, "]", sep = "")]
  alpha.13 <- params[,paste("alpha.13[", s, "]", sep = "")]
  alpha.22 <- params[,paste("alpha.22[", s, "]", sep = "")]
  alpha.33 <- params[,paste("alpha.33[", s, "]", sep = "")]
  alpha.32 <- params[,paste("alpha.32[", s, "]", sep = "")]
  
  # select appropriate parameters
  if("parameter" %in% type){
    phi.l.mu <- params[draw,"phi.l.mu"]
    phi.n.mu <- params[draw,"phi.n.mu"]
    phi.a.mu <- params[draw,"phi.a.mu"]
    grow.ln.mu <- params[draw,"grow.ln.mu"]
    grow.na.mu <- params[draw,"grow.na.mu"]
    repro.mu <- params[draw,"repro.mu"]
    beta.11 <- params[draw, "beta.11"]
    beta.22 <- params[draw, "beta.22"]
    beta.33 <- params[draw, "beta.33"]
    beta.R.11 <- molt.nymph[draw,"beta"]
    intercept.R.11 <- molt.nymph[draw,"lambda.mu"]
    alpha.11 <- alpha.11[draw]
    alpha.13 <- alpha.13[draw]
    alpha.22 <- alpha.22[draw]
    alpha.33 <- alpha.33[draw]
    alpha.32 <- alpha.32[draw]
  } else {
    phi.l.mu <- rep(param.mean["phi.l.mu"], Nmc)
    phi.n.mu <- rep(param.mean["phi.n.mu"], Nmc)
    phi.a.mu <- rep(param.mean["phi.a.mu"], Nmc)
    grow.ln.mu <- rep(param.mean["grow.ln.mu"], Nmc)
    grow.na.mu <- rep(param.mean["grow.na.mu"], Nmc)
    repro.mu <- rep(param.mean["repro.mu"], Nmc)
    beta.11 <- rep(param.mean["beta.11"], Nmc)
    beta.22 <- rep(param.mean["beta.22"], Nmc)
    beta.33 <- rep(param.mean["beta.33"], Nmc)
    beta.R.11 <- rep(molt.nymph.surv.mean["beta"], Nmc)
    intercept.R.11 <- rep(molt.nymph.surv.mean["lambda.mu"], Nmc)
    alpha.11 <- rep(mean(alpha.11), Nmc)
    alpha.13 <- rep(mean(alpha.13), Nmc)
    alpha.22 <- rep(mean(alpha.22), Nmc)
    alpha.33 <- rep(mean(alpha.33), Nmc)
    alpha.32 <- rep(mean(alpha.32), Nmc)
  }
  
  # select appropraite covariance matrix
  if("process" %in% type){
    SIGMA <- array(NA, dim = c(3,3,Nmc))
    SIGMA[1,1,] <- params[draw,"SIGMA[1,1]"]
    SIGMA[1,2,] <- params[draw,"SIGMA[1,2]"]
    SIGMA[1,3,] <- params[draw,"SIGMA[1,3]"]
    SIGMA[2,1,] <- params[draw,"SIGMA[2,1]"]
    SIGMA[2,2,] <- params[draw,"SIGMA[2,2]"]
    SIGMA[2,3,] <- params[draw,"SIGMA[2,3]"]
    SIGMA[3,1,] <- params[draw,"SIGMA[3,1]"]
    SIGMA[3,2,] <- params[draw,"SIGMA[3,2]"]
    SIGMA[3,3,] <- params[draw,"SIGMA[3,3]"]
    
    # convert from precision to standard dev
    for(i in 1:Nmc){
      SIGMA[,,i] <- solve(SIGMA[,,i])
    }
  } else if("process mean" %in% type){
    SIGMA <- array(NA, dim = c(3,3,Nmc))
    SIGMA[1,1,] <- rep(param.mean["SIGMA[1,1]"], Nmc)
    SIGMA[1,2,] <- rep(param.mean["SIGMA[1,2]"], Nmc)
    SIGMA[1,3,] <- rep(param.mean["SIGMA[1,3]"], Nmc)
    SIGMA[2,1,] <- rep(param.mean["SIGMA[2,1]"], Nmc)
    SIGMA[2,2,] <- rep(param.mean["SIGMA[2,2]"], Nmc)
    SIGMA[2,3,] <- rep(param.mean["SIGMA[2,3]"], Nmc)
    SIGMA[3,1,] <- rep(param.mean["SIGMA[3,1]"], Nmc)
    SIGMA[3,2,] <- rep(param.mean["SIGMA[3,2]"], Nmc)
    SIGMA[3,3,] <- rep(param.mean["SIGMA[3,3]"], Nmc)
    
    # convert from precision to standard dev
    for(i in 1:Nmc){
      SIGMA[,,i] <- solve(SIGMA[,,i])
    }
  } else {
    SIGMA <- array(0, dim = c(3,3,Nmc))
  }
  
  # random error
  if("random effect" %in% type){
    tau.11 <- 1/sqrt(params[draw, "tau.11"])
    tau.13 <- 1/sqrt(params[draw, "tau.13"])
    tau.22 <- 1/sqrt(params[draw, "tau.22"])
    tau.33 <- 1/sqrt(params[draw, "tau.33"])
    tau.32 <- 1/sqrt(params[draw, "tau.32"])
    alpha.11 <- rnorm(Nmc, 0, tau.11)
    alpha.13 <- rnorm(Nmc, 0, tau.13)
    alpha.22 <- rnorm(Nmc, 0, tau.22)
    alpha.33 <- rnorm(Nmc, 0, tau.33)
    alpha.32 <- rnorm(Nmc, 0, tau.32)
  }
  
  # run mcmc sampling
  for(m in 1:Nmc){
    
    # draw transition matrix
    A[1,1,] <- inv.logit(phi.l.mu[m] + beta.11[m]*temp + alpha.11[m])
    A[2,2,] <- inv.logit(phi.n.mu[m] + beta.22[m]*temp + alpha.22[m])
    A[3,3,] <- inv.logit(phi.a.mu[m] + beta.33[m]*temp + alpha.33[m])
    
    for(r in 1:N_days){
      P[1,1,r] <- rbinom(1,1,inv.logit(0.07*cum.gdd[r]))       # egg to questing larva
      Q[2,2,r] <- rbinom(1,1,inv.logit(0.0004*cum.gdd[r]))     # molting adult to questing adult
      S[1,1,r] <- rbinom(1,1,inv.logit(-0.009*cum.gdd[r]))     # molting nymph to questing nymph
    }
    # spring-to-summer matrix
    P[2,2,] <- inv.logit(grow.na.mu[m] + alpha.32[m])  # questing nymph to molting adult
    
    # summer-to-fall matrix
    Q[1,1,] <- inv.logit(grow.ln.mu[m])          # questing larva to molting nymph
    
    # fall-to-winter matrix
    R[1,1,] <- inv.logit(intercept.R.11[m] + beta.R.11[m]*temp) # molting nymph survival
    R[2,2,] <- rep(q.adult_d.adult, N_days)                     # questing adult to dormant adult transition
    
    # winter-to-spring matrix
    
    S[2,2,] <- exp(repro.mu[m] + alpha.13[m])      # dormant adult to egg (reproduction)
    
    e <- IC[1]       # eggs
    l <- IC[1]       # questing larvae
    m.n <- IC[2]     # molting nymphs
    n <- IC[2]       # questing nymphs
    m.a <- IC[2]     # molting adults
    a <- IC[3]       # questing adults
    d.a <- IC[3]     # dormant adutls

    ## aggrigate transition matricies
    for(t in 1:(N_est)){   # loop over the number of sampling days
      if(t == 1){
        A.TRANS <- A[,,1] %*% A[,,2]
        for(d in 3:(df[1]-2)){
          A.TRANS <- A.TRANS %*% A[,,d]
        }
      } else {
        for(d in 1:df[t]){
          if(d == 1){
            A.TRANS <- A[,,dt.index[t-1]] %*% A[,,dt.index[t-1]+1]
          } else {
            A.TRANS <- A.TRANS %*% A[,,dt.index[t-1]+d-1]
          }
        }
      }
      
      # predict questing ticks
      quest.process <- A.TRANS %*% c(l,n,a)
      
      # questing process error
      Ex.quest <- rmvnorm(1,quest.process,SIGMA[,,m]) 
      
      # predict seasonal transitions
      spring_summer <- P[,,dt.index[t]] %*% c(e, n)
      summer_fall <- Q[,,dt.index[t]] %*% c(l, m.a)
      fall_winter <- R[,,dt.index[t]] %*% c(m.n, a)
      winter_spring <- S[,,dt.index[t]] %*% c(m.n, d.a)
      
      # extract values
      e <- winter_spring[2]                           # eggs
      l <- max(Ex.quest[1], 0) + spring_summer[1]     # questing larvae
      m.n <- summer_fall[1] + fall_winter[1]          # molting nymphs
      n <- max(Ex.quest[2], 0) + winter_spring[1]     # questing nymphs
      m.a <- spring_summer[2]                         # molting adults
      a <- max(Ex.quest[3], 0) + summer_fall[2]       # questing adults
      d.a <- fall_winter[2]                           # dormant adutls

      # store
      pred[,t,m] <- c(e,l,m.n,n,m.a,a,d.a) 
    }
  }
  return(pred)
} # close function

# egg2larva <- 0.1
# nymph2adult <- 0.05
# larva2nymph <- 0.05
# molt.adult2adult <- 0.3
q.adult_d.adult <- 1
# molt.nymph.surv <- 0.8
# molt.nymph2nymph <- 0.98

Nmc <- 500 
draw <- sample.int(nrow(params.mat), Nmc, replace = TRUE)
N_est <- 20
df <- 21

pred <- forecast_with_molting_temp(type = c("process"),
                                     site = "Henry Control",
                                     params = params.mat,
                                     ic = ic.mat,
                                     N_est = N_est,
                                     df = df,
                                     Nmc = Nmc,
                                     draw = draw,
                                     q.adult_d.adult,molt.nymph)

egg.pred <- apply(pred[1,,], 1, quantile, c(0.025,0.5,0.975))
larv.pred <- apply(pred[2,,], 1, quantile, c(0.025,0.5,0.975))
molt.nymph.pred <- apply(pred[3,,], 1, quantile, c(0.025,0.5,0.975))
nymph.pred <- apply(pred[4,,], 1, quantile, c(0.025,0.5,0.975))
molt.adult.pred <- apply(pred[5,,], 1, quantile, c(0.025,0.5,0.975))
adult.pred <- apply(pred[6,,], 1, quantile, c(0.025,0.5,0.975))
dorm.adult.pred <- apply(pred[7,,], 1, quantile, c(0.025,0.5,0.975))
time <- 1:ncol(larv.pred)


s <- 2

last <- get_last_day()
last.day <- as.Date(last[s,"last.day"])

date.seq <- seq.Date(last.day, by = df, length.out = N_est)
ym.seq <- format(date.seq, "%Y-%m")
at <- seq(1, N_est, length.out = 7)

plot(time,larv.pred[2,],pch="",main="Larvae", xaxt="n")
ciEnvelope(time,larv.pred[1,],larv.pred[3,],col="lightblue")
lines(time,larv.pred[2,])
axis(1, at = at, labels = ym.seq[at])

plot(time,nymph.pred[2,],pch="",main="Nymph", xaxt="n")
ciEnvelope(time,nymph.pred[1,],nymph.pred[3,],col="lightblue")
lines(time,nymph.pred[2,])
axis(1, at = at, labels = ym.seq[at])

plot(time,adult.pred[2,],pch="",main="Adult", xaxt="n")
ciEnvelope(time,adult.pred[1,],adult.pred[3,],col="lightblue")
lines(time,adult.pred[2,])
axis(1, at = at, labels = ym.seq[at])

# sum(which(nymph.pred==adult.pred))

stage <- c("Egg","Larvae","Molting Nymph","Nymph","Molting Adult","Adult","Dormant Adult")
par(mfrow = c(1,1))
plot(time, egg.pred[2,], type = "l", ylim = c(0,1000), xaxt="n")
lines(time, larv.pred[2,], col = 2)
lines(time, molt.nymph.pred[2,], col = 3)
lines(time, nymph.pred[2,], col = 4)
lines(time, molt.adult.pred[2,], col = 5)
lines(time, adult.pred[2,], col = 6)
lines(time, dorm.adult.pred[2,], col = 7)
axis(1, at = at, labels = ym.seq[at])
legend("topright",legend = stage,lty = 1,col=1:7)

