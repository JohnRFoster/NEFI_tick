library(ecoforecastR)

source("Models/predict_state_one_gdd_null.R")
source("Models/forecast_state_gdd_null.R")
source("Functions/combine_chains.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/life_stage_ci.R")
source("Functions/site_data_met.R")

site <- "Green Control"

data <- cary_ticks_met_JAGS()
data.site <- site_data_met(site,NULL,data)
path <- "../FinalOut/Independent_Fits/GDDThreshold/LowThreshOnly/GDDSwitch_GreenControl"

out <- combine_chains(path,1:5)

params <- as.matrix(out$params)
ic <- as.matrix(out$predict)
# ic <- out$predict.m

N_est <- 30
df <- 21

pred.D <- predict_state_one_gdd_null("deterministic",
                                     "low",
                                     site,
                                     params,
                                     ic,
                                     data,
                                     1,1)

forecast.D <- forecast_state_gdd_null("deterministic",
                                     "low",
                                     site,
                                     params,
                                     ic,
                                     N_est,
                                     df,
                                     1,1)


type="deterministic"
thresh="low"
site=site
params=params
ic=ic
data=data
Nmc=1
draw=1

# plot(1:72,pred.D[1,,1],type="l",ylim=c(0,3000))
# points(1:72,data.site$y[1,])
# 
# plot(1:72,pred.D[2,,1],type="l",ylim=c(0,100))
# points(1:72,data.site$y[2,])
# 
# plot(1:72,pred.D[3,,1],type="l",ylim=c(0,40))
# points(1:72,data.site$y[3,])

days <- 72+N_est

param.col <- "gray30"
param.ic.col <- "honeydew4"
param.ic.drive.col <- "honeydew3"
param.ic.drive.process.col <- "honeydew2"

Nmc <- 500 
draw <- sample.int(nrow(params), Nmc)


### Parameter Unertainty ###

predict.param <- predict_state_one_gdd_null("parameter",
                                     "low",
                                     site,
                                     params,
                                     ic,
                                     data,
                                     Nmc,draw)

forecast.param <- forecast_state_gdd_null("parameter",
                                      "low",
                                      site,
                                      params,
                                      ic,
                                      N_est,
                                      df,
                                      Nmc,draw)

predict.95ci.P <- life.stage.ci(predict.param,"predict")
forecast.95ci.P <- life.stage.ci(forecast.param,"predict")

time.1 <- 1:ncol(predict.95ci.P[[1]])
time.2 <- (ncol(predict.95ci.P[[1]])+1):(days-2)
time.all <- 1:(days-1)

plot(time.all,pch="",ylim=c(0,2000),main="Larvae")
ciEnvelope(time.1,predict.95ci.P[[1]][1,],predict.95ci.P[[1]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[1]][1,],forecast.95ci.P[[1]][3,],col=param.col)
lines(time.1,predict.95ci.P[[1]][2,])
lines(time.2,forecast.95ci.P[[1]][2,])

plot(time.all,pch="",ylim=c(0,100),main="Nymph")
ciEnvelope(time.1,predict.95ci.P[[2]][1,],predict.95ci.P[[2]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[2]][1,],forecast.95ci.P[[2]][3,],col=param.col)
lines(time.1,predict.95ci.P[[2]][2,])
lines(time.2,forecast.95ci.P[[2]][2,])

plot(time.all,pch="",ylim=c(0,40),main="Adult")
ciEnvelope(time.1,predict.95ci.P[[3]][1,],predict.95ci.P[[3]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[3]][1,],forecast.95ci.P[[3]][3,],col=param.col)
lines(time.1,predict.95ci.P[[3]][2,])
lines(time.2,forecast.95ci.P[[3]][2,])


### Parameter + Initial Condition Unertainty ###
type <- c("parameter","ic")

predict.ic.param <- predict_state_one_gdd_null(type,
                                         "low",
                                         site,
                                         params,
                                         ic,
                                         data,
                                         Nmc,draw)

forecast.ic.param <- forecast_state_gdd_null(type,
                                       "low",
                                       site,
                                       params,
                                       ic,
                                       N_est,
                                       df,
                                       Nmc,draw)

predict.95ci.P.IC <- life.stage.ci(predict.ic.param,"predict")
forecast.95ci.P.IC <- life.stage.ci(forecast.ic.param,"predict")

plot(time.all,pch="",ylim=c(0,2000),main="Larvae")
ciEnvelope(time.1,predict.95ci.P.IC[[1]][1,],predict.95ci.P.IC[[1]][3,],col=param.ic.col)
ciEnvelope(time.2,forecast.95ci.P.IC[[1]][1,],forecast.95ci.P.IC[[1]][3,],col=param.ic.col)
ciEnvelope(time.1,predict.95ci.P[[1]][1,],predict.95ci.P[[1]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[1]][1,],forecast.95ci.P[[1]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[1]][2,])
lines(time.2,forecast.95ci.P.IC[[1]][2,])

plot(time.all,pch="",ylim=c(0,100),main="Nymph")
ciEnvelope(time.1,predict.95ci.P.IC[[2]][1,],predict.95ci.P.IC[[2]][3,],col=param.ic.col)
ciEnvelope(time.2,forecast.95ci.P.IC[[2]][1,],forecast.95ci.P.IC[[2]][3,],col=param.ic.col)
ciEnvelope(time.1,predict.95ci.P[[2]][1,],predict.95ci.P[[2]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[2]][1,],forecast.95ci.P[[2]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[2]][2,])
lines(time.2,forecast.95ci.P[[2]][2,])

plot(time.all,pch="",ylim=c(0,40),main="Adult")
ciEnvelope(time.1,predict.95ci.P.IC[[3]][1,],predict.95ci.P.IC[[3]][3,],col=param.ic.col)
ciEnvelope(time.2,forecast.95ci.P.IC[[3]][1,],forecast.95ci.P.IC[[3]][3,],col=param.ic.col)
ciEnvelope(time.1,predict.95ci.P[[3]][1,],predict.95ci.P[[3]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[3]][1,],forecast.95ci.P[[3]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[3]][2,])
lines(time.2,forecast.95ci.P.IC[[3]][2,])


### Process + Parameter + Initial Condition Unertainty ###
type <- c("process","parameter","ic")

predict.ic.param.proc <- predict_state_one_gdd_null(type,
                                         "low",
                                         site,
                                         params,
                                         ic,
                                         data,
                                         Nmc,draw)

forecast.ic.param.proc <- forecast_state_gdd_null(type,
                                       "low",
                                       site,
                                       params,
                                       ic,
                                       N_est,
                                       df,
                                       Nmc,draw)

predict.95ci.Proc.P.IC <- life.stage.ci(predict.ic.param.proc,"predict")
forecast.95ci.Proc.P.IC <- life.stage.ci(forecast.ic.param.proc,"predict")

plot(time.all,pch="",ylim=c(0,2000),main="Larvae")
ciEnvelope(time.1,predict.95ci.Proc.P.IC[[1]][1,],predict.95ci.Proc.P.IC[[1]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.2,forecast.95ci.Proc.P.IC[[1]][1,],forecast.95ci.Proc.P.IC[[1]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.1,predict.95ci.P.IC[[1]][1,],predict.95ci.P.IC[[1]][3,],col=param.ic.col)
ciEnvelope(time.2,forecast.95ci.P.IC[[1]][1,],forecast.95ci.P.IC[[1]][3,],col=param.ic.col)
ciEnvelope(time.1,predict.95ci.P[[1]][1,],predict.95ci.P[[1]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[1]][1,],forecast.95ci.P[[1]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[1]][2,])
lines(time.2,forecast.95ci.P.IC[[1]][2,])

plot(time.all,pch="",ylim=c(0,100),main="Nymph")
ciEnvelope(time.1,predict.95ci.Proc.P.IC[[2]][1,],predict.95ci.Proc.P.IC[[2]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.2,forecast.95ci.Proc.P.IC[[2]][1,],forecast.95ci.Proc.P.IC[[2]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.1,predict.95ci.P.IC[[2]][1,],predict.95ci.P.IC[[2]][3,],col=param.ic.col)
ciEnvelope(time.2,forecast.95ci.P.IC[[2]][1,],forecast.95ci.P.IC[[2]][3,],col=param.ic.col)
ciEnvelope(time.1,predict.95ci.P[[2]][1,],predict.95ci.P[[2]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[2]][1,],forecast.95ci.P[[2]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[2]][2,])
lines(time.2,forecast.95ci.P[[2]][2,])

plot(time.all,pch="",ylim=c(0,40),main="Adult")
ciEnvelope(time.1,predict.95ci.Proc.P.IC[[3]][1,],predict.95ci.Proc.P.IC[[3]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.2,forecast.95ci.Proc.P.IC[[3]][1,],forecast.95ci.Proc.P.IC[[3]][3,],col=param.ic.drive.process.col)
ciEnvelope(time.1,predict.95ci.P.IC[[3]][1,],predict.95ci.P.IC[[3]][3,],col=param.ic.col)
ciEnvelope(time.2,forecast.95ci.P.IC[[3]][1,],forecast.95ci.P.IC[[3]][3,],col=param.ic.col)
ciEnvelope(time.1,predict.95ci.P[[3]][1,],predict.95ci.P[[3]][3,],col=param.col)
ciEnvelope(time.2,forecast.95ci.P[[3]][1,],forecast.95ci.P[[3]][3,],col=param.col)
lines(time.1,predict.95ci.P.IC[[3]][2,])
lines(time.2,forecast.95ci.P.IC[[3]][2,])

