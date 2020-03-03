library(mvtnorm)
library(ecoforecastR)
library(boot)

source("Functions/Future_Met.R")
source("Functions/life_stage_ci.R")

source("Models/predict_state_one_gdd_K_set_independent.R")
source("Models/forecast_state_one_gdd_K_set_independent.R")
source("Functions/ammend_chains.R")
source("Functions/cary_tick_met_JAGS.R")


# dir <- "../FinalOut/HB_Partial_GDD/Temp/LowLoop/GDDSwitch_Low_Temp"
# dir <- "../FinalOut/HB_Partial_GDD/K_estimate/WindowLoop_4alpha/GDDSwitch_K_Window4alpha"
# 
# out <- ammend_chains(dir, 5, 24, save = TRUE)

load("../FinalOut/Independent_Fits/GDDThreshold/K_all/K_set_multivariate/Green_normTransition/Combined_thinMat_K_set_GreenControl.RData")

data <- cary_ticks_met_JAGS()

params <- as.matrix(out.test$params)
ic <- as.matrix(out.test$predict)
n.iter <- nrow(params.mat)

Nmc <- 250 
draw <- sample.int(nrow(params.mat), Nmc)
N_est <- 40
df <- 16
site <- "Green Control"

# pred.D <- predict_state_one_gdd_k_window("deterministic",
#                                      site,
#                                      params,
#                                      ic,
#                                      data,
#                                      1,1)
# 
# forecast.D <- forecast_state_one_gdd_k_window("deterministic",
#                                       site,
#                                       params,
#                                       ic,
#                                       N_est,
#                                       df,
#                                       1,1)
# 
days <- 72+N_est

param.col <- "gray30"
param.ic.col <- "honeydew4"
param.ic.drive.col <- "honeydew3"
param.ic.drive.process.col <- "honeydew2"
# 
# ### Parameter Unertainty ###
# 
# predict.param <- predict_state_one_gdd_k_window("parameter",
#                                             site,
#                                             params,
#                                             ic,
#                                             data,
#                                             Nmc,draw)
# 
# forecast.param <- forecast_state_one_gdd_k_window("parameter",
#                                           site,
#                                           params,
#                                           ic,
#                                           N_est,
#                                           df,
#                                           Nmc,draw)
# 
# predict.95ci.P <- life.stage.ci(predict.param,"predict")
# forecast.95ci.P <- life.stage.ci(forecast.param,"predict")
# 
# time.1 <- 1:ncol(predict.95ci.P[[1]])
# time.2 <- (ncol(predict.95ci.P[[1]])+1):(days-1)
# time.all <- 1:(days-1)
# 
# plot(time.all,pch="",ylim=c(0,2000),main="Larvae")
# ciEnvelope(time.1,predict.95ci.P[[1]][1,],predict.95ci.P[[1]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[1]][1,],forecast.95ci.P[[1]][3,],col=param.col)
# lines(time.1,predict.95ci.P[[1]][2,])
# lines(time.2,forecast.95ci.P[[1]][2,])
# 
# plot(time.all,pch="",ylim=c(0,100),main="Nymph")
# ciEnvelope(time.1,predict.95ci.P[[2]][1,],predict.95ci.P[[2]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[2]][1,],forecast.95ci.P[[2]][3,],col=param.col)
# lines(time.1,predict.95ci.P[[2]][2,])
# lines(time.2,forecast.95ci.P[[2]][2,])
# 
# plot(time.all,pch="",ylim=c(0,40),main="Adult")
# ciEnvelope(time.1,predict.95ci.P[[3]][1,],predict.95ci.P[[3]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[3]][1,],forecast.95ci.P[[3]][3,],col=param.col)
# lines(time.1,predict.95ci.P[[3]][2,])
# lines(time.2,forecast.95ci.P[[3]][2,])
# 
# 
# ### Parameter + Initial Condition Unertainty ###
# type <- c("parameter","ic")
# 
# predict.ic.param <- predict_state_one_gdd_k_window(type,
#                                                site,
#                                                params,
#                                                ic,
#                                                data,
#                                                Nmc,draw)
# 
# forecast.ic.param <- forecast_state_one_gdd_k_window(type,
#                                              site,
#                                              params,
#                                              ic,
#                                              N_est,
#                                              df,
#                                              Nmc,draw)
# 
# predict.95ci.P.IC <- life.stage.ci(predict.ic.param,"predict")
# forecast.95ci.P.IC <- life.stage.ci(forecast.ic.param,"predict")
# 
# plot(time.all,pch="",ylim=c(0,2000),main="Larvae")
# ciEnvelope(time.1,predict.95ci.P.IC[[1]][1,],predict.95ci.P.IC[[1]][3,],col=param.ic.col)
# ciEnvelope(time.2,forecast.95ci.P.IC[[1]][1,],forecast.95ci.P.IC[[1]][3,],col=param.ic.col)
# ciEnvelope(time.1,predict.95ci.P[[1]][1,],predict.95ci.P[[1]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[1]][1,],forecast.95ci.P[[1]][3,],col=param.col)
# lines(time.1,predict.95ci.P.IC[[1]][2,])
# lines(time.2,forecast.95ci.P.IC[[1]][2,])
# 
# plot(time.all,pch="",ylim=c(0,100),main="Nymph")
# ciEnvelope(time.1,predict.95ci.P.IC[[2]][1,],predict.95ci.P.IC[[2]][3,],col=param.ic.col)
# ciEnvelope(time.2,forecast.95ci.P.IC[[2]][1,],forecast.95ci.P.IC[[2]][3,],col=param.ic.col)
# ciEnvelope(time.1,predict.95ci.P[[2]][1,],predict.95ci.P[[2]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[2]][1,],forecast.95ci.P[[2]][3,],col=param.col)
# lines(time.1,predict.95ci.P.IC[[2]][2,])
# lines(time.2,forecast.95ci.P[[2]][2,])
# 
# plot(time.all,pch="",ylim=c(0,40),main="Adult")
# ciEnvelope(time.1,predict.95ci.P.IC[[3]][1,],predict.95ci.P.IC[[3]][3,],col=param.ic.col)
# ciEnvelope(time.2,forecast.95ci.P.IC[[3]][1,],forecast.95ci.P.IC[[3]][3,],col=param.ic.col)
# ciEnvelope(time.1,predict.95ci.P[[3]][1,],predict.95ci.P[[3]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[3]][1,],forecast.95ci.P[[3]][3,],col=param.col)
# lines(time.1,predict.95ci.P.IC[[3]][2,])
# lines(time.2,forecast.95ci.P.IC[[3]][2,])


### Process + Parameter + Initial Condition Unertainty ###
type <- c("process","parameter","ic")

predict.ic.param.proc <- predict_state_one_gdd_K_set_independent(type,
                                                    site,
                                                    params.mat,
                                                    predict.mat,
                                                    data,
                                                    Nmc,draw)

forecast.ic.param.proc <- forecast_state_one_gdd_K_set_independent(type,
                                                  site,
                                                  params.mat,
                                                  predict.mat,
                                                  N_est,
                                                  df,
                                                  Nmc,draw)

predict.95ci.Proc.P.IC <- life.stage.ci(predict.ic.param.proc,"predict")
forecast.95ci.Proc.P.IC <- life.stage.ci(forecast.ic.param.proc,"predict")

time.1 <- 1:ncol(predict.95ci.Proc.P.IC[[1]])
time.2 <- time.1[length(time.1)] + 1:ncol(forecast.95ci.Proc.P.IC[[1]])
time.all <- 1:time.2[length(time.2)]

plot(time.all,pch="",ylim=c(0,2000),main="Larvae")
ciEnvelope(time.1,predict.95ci.Proc.P.IC[[1]][1,],predict.95ci.Proc.P.IC[[1]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.2,forecast.95ci.Proc.P.IC[[1]][1,],forecast.95ci.Proc.P.IC[[1]][3,],col=param.ic.drive.process.col)
# ciEnvelope(time.1,predict.95ci.P.IC[[1]][1,],predict.95ci.P.IC[[1]][3,],col=param.ic.col)
# ciEnvelope(time.2,forecast.95ci.P.IC[[1]][1,],forecast.95ci.P.IC[[1]][3,],col=param.ic.col)
# ciEnvelope(time.1,predict.95ci.P[[1]][1,],predict.95ci.P[[1]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[1]][1,],forecast.95ci.P[[1]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[1]][2,])
lines(time.2,forecast.95ci.P.IC[[1]][2,])
lines(time.2,forecast.95ci.Proc.P.IC[[1]][2,], col = 2)

plot(time.all,pch="",ylim=c(0,250),main="Nymph")
ciEnvelope(time.1,predict.95ci.Proc.P.IC[[2]][1,],predict.95ci.Proc.P.IC[[2]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.2,forecast.95ci.Proc.P.IC[[2]][1,],forecast.95ci.Proc.P.IC[[2]][3,],col=param.ic.drive.process.col)
# ciEnvelope(time.1,predict.95ci.P.IC[[2]][1,],predict.95ci.P.IC[[2]][3,],col=param.ic.col)
# ciEnvelope(time.2,forecast.95ci.P.IC[[2]][1,],forecast.95ci.P.IC[[2]][3,],col=param.ic.col)
# ciEnvelope(time.1,predict.95ci.P[[2]][1,],predict.95ci.P[[2]][3,],col=param.col)
# ciEnvelope(time.2,forecast.95ci.P[[2]][1,],forecast.95ci.P[[2]][3,],col=param.col)
lines(time.1,predict.95ci.Proc.P.IC[[2]][2,])
# lines(time.2,forecast.95ci.P[[2]][2,])
lines(time.2,forecast.95ci.Proc.P.IC[[2]][2,], col = 2)

plot(time.all,pch="",ylim=c(0,40),main="Adult")
ciEnvelope(time.1,predict.95ci.Proc.P.IC[[3]][1,],predict.95ci.Proc.P.IC[[3]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.2,forecast.95ci.Proc.P.IC[[3]][1,],forecast.95ci.Proc.P.IC[[3]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.1,predict.95ci.P.IC[[3]][1,],predict.95ci.P.IC[[3]][3,],col=param.ic.col)
ciEnvelope(time.2,forecast.95ci.P.IC[[3]][1,],forecast.95ci.P.IC[[3]][3,],col=param.ic.col)
ciEnvelope(time.1,predict.95ci.P[[3]][1,],predict.95ci.P[[3]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[3]][1,],forecast.95ci.P[[3]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[3]][2,])
lines(time.2,forecast.95ci.P.IC[[3]][2,])
lines(time.2,forecast.95ci.Proc.P.IC[[3]][2,], col = 2)