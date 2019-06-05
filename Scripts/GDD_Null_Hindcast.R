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

ic.col <- "gray30"
param.col <- "honeydew4"
drive.col <- "honeydew3"
process.col <- "honeydew2"

Nmc <- 500 
draw <- sample.int(nrow(params), Nmc)


### Initial Condition Unertainty ###

predict.IC <- predict_state_one_gdd_null("ic",
                                     "low",
                                     site,
                                     params,
                                     ic,
                                     data,
                                     Nmc,draw)

forecast.IC <- forecast_state_gdd_null("ic",
                                      "low",
                                      site,
                                      params,
                                      ic,
                                      N_est,
                                      df,
                                      Nmc,draw)

predict.95ci.IC <- life.stage.ci(predict.IC,"predict")
forecast.95ci.IC <- life.stage.ci(forecast.IC,"predict")

time.1 <- 1:ncol(predict.95ci.IC[[1]])
time.2 <- (ncol(predict.95ci.IC[[1]])+1):(days-2)
time.all <- 1:(days-1)

plot(time.all,pch="",ylim=c(0,2000),main="Larvae")
ciEnvelope(time.1,predict.95ci.IC[[1]][1,],predict.95ci.IC[[1]][3,],col=ic.col)
ciEnvelope(time.2,forecast.95ci.IC[[1]][1,],forecast.95ci.IC[[1]][3,],col=ic.col)
lines(time.1,predict.95ci.IC[[1]][2,])
lines(time.2,forecast.95ci.IC[[1]][2,])

plot(time.all,pch="",ylim=c(0,100),main="Nymph")
ciEnvelope(time.1,predict.95ci.IC[[2]][1,],predict.95ci.IC[[2]][3,],col=ic.col)
ciEnvelope(time.2,forecast.95ci.IC[[2]][1,],forecast.95ci.IC[[2]][3,],col=ic.col)
lines(time.1,predict.95ci.IC[[2]][2,])
lines(time.2,forecast.95ci.IC[[2]][2,])

plot(time.all,pch="",ylim=c(0,40),main="Adult")
ciEnvelope(time.1,predict.95ci.IC[[3]][1,],predict.95ci.IC[[3]][3,],col=ic.col)
ciEnvelope(time.2,forecast.95ci.IC[[3]][1,],forecast.95ci.IC[[3]][3,],col=ic.col)
lines(time.1,predict.95ci.IC[[3]][2,])
lines(time.2,forecast.95ci.IC[[3]][2,])




