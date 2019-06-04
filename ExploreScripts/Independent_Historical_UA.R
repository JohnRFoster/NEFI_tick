
library(knitr)
library(ecoforecastR)
library(runjags)
library(fosteR)
library(mvtnorm)
library(ggplot2)
library(grid)
library(gridExtra)
library(boot)

source("Functions/predict_state_one_null.R")
source("Functions/predict_state_one_precip.R")
source("Functions/gg_data_ci_partition.R")
source("Functions/model_performance_metrics.R")

life.stage.ci <- function(pred, quants = c(0.025, 0.5, 0.975)){
  time <- 1:(length(pred[1,,1])-1)
  ci <- list()
  larv <- pred[1,,]
  ci[[1]] <- apply(larv[time,], 1, quantile, quants)
  
  nymph <- pred[2,,]
  ci[[2]] <- apply(nymph[time,], 1, quantile, quants)
  
  adult <- pred[3,,]
  ci[[3]] <- apply(adult[time,], 1, quantile, quants)
  
  names(ci) <- c("Larvae","Nymph","Adult")
  return(ci)
}

load("/projectnb/dietzelab/fosterj/FinalOut/HB_Partial_Precip/Out_Precip_Partial_HB.RData")
summary(params)[[1]]

params.mat <- as.matrix(params)
ic.mat <- as.matrix(predict)
# params <- window(params, start = 3000)
# predict <- window(predict, start = 3000)
n.iter <- nrow(params.mat)

Nmc <- 500 
draw <- sample.int(nrow(params.mat), Nmc)

sites <- c("Green Control","Henry Control","Tea Control")
data <- cary_ticks_met_JAGS(sites)
N_days <- data$N_days
precip <- data$met[,3,]
nymph.sum <- larv.sum <- data.frame()

for(i in 1:3){
  x <- precip[1:N_days[i],i]
  l.sum <- matrix(NA, Nmc, length(x))
  n.sum <- matrix(NA, Nmc, length(x))
  for(t in 1:N_days[i]){
    l.sum[,t] <- inv.logit(params.mat[draw,"phi.l.mu"] + params.mat[draw,"beta.11"]*x[t]) + inv.logit(params.mat[draw,"grow.ln.mu"])
    n.sum[,t] <- inv.logit(params.mat[draw,"phi.n.mu"] + params.mat[draw,"beta.22"]*x[t] + params.mat[draw,paste("alpha.22[", i, "]", sep = "")]) + 
      inv.logit(params.mat[draw,"grow.na.mu"] + params.mat[draw,paste("alpha.32[", i, "]", sep = "")])
  }
  larv.sum[i,1] <- length(which(l.sum > 1)) / (Nmc * N_days[i]) * 100
  larv.sum[i,2] <- max(l.sum)
  
  nymph.sum[i,1] <- length(which(n.sum > 1)) / (Nmc * N_days[i]) * 100
  nymph.sum[i,2] <- max(n.sum)
}

larv.sum
nymph.sum

n.sum <- apply(n.sum, 2, quantile, c(0.025,0.5,0.95))

mean(n.sum[2,])


pred.param <- list()
pred.param$green <- predict_state_one_precip(type = "parameter",
                                           site = "Green Control",
                                           params = params.mat,
                                           ic = ic.mat,
                                           data = data,
                                           Nmc = Nmc,
                                           draw = draw) 
pred.param$henry <- predict_state_one_precip(type = "parameter",
                                           site = "Henry Control",
                                           params = params.mat,
                                           ic = ic.mat,
                                           data = data,
                                           Nmc = Nmc,
                                           draw = draw) 
pred.param$tea <- predict_state_one_precip(type = "parameter",
                                         site = "Tea Control",
                                         params = params.mat,
                                         ic = ic.mat,
                                         data = data,
                                         Nmc = Nmc,
                                         draw = draw) 
conf.int.param <- list()
conf.int.param$green <- life.stage.ci(pred.param$green)
conf.int.param$henry <- life.stage.ci(pred.param$henry)
conf.int.param$tea <- life.stage.ci(pred.param$tea)

# param + IC
pred.ic.param <- list()
type <- c("parameter", "ic")
pred.ic.param$green <- predict_state_one_precip(type = type,
                                              site = "Green Control",
                                              params = params.mat,
                                              ic = ic.mat,
                                              data = data,
                                              Nmc = Nmc,
                                              draw = draw) 
pred.ic.param$henry <- predict_state_one_precip(type = type,
                                              site = "Henry Control",
                                              params = params.mat,
                                              ic = ic.mat,
                                              data = data,
                                              Nmc = Nmc,
                                              draw = draw) 
pred.ic.param$tea <- predict_state_one_precip(type = type,
                                            site = "Tea Control",
                                            params = params.mat,
                                            ic = ic.mat,
                                            data = data,
                                            Nmc = Nmc,
                                            draw = draw)

conf.int.ic.param <- list()
conf.int.ic.param$green <- life.stage.ci(pred.ic.param$green)
conf.int.ic.param$henry <- life.stage.ci(pred.ic.param$henry)
conf.int.ic.param$tea <- life.stage.ci(pred.ic.param$tea)

# param + IC + Process
pred.ic.param.process <- list()
type <- c("parameter", "ic", "process")
pred.ic.param.process$green <- predict_state_one_precip(type = type,
                                                      site = "Green Control",
                                                      params = params.mat,
                                                      ic = ic.mat,
                                                      data = data,
                                                      Nmc = Nmc,
                                                      draw = draw) 
pred.ic.param.process$henry <- predict_state_one_precip(type = type,
                                                      site = "Henry Control",
                                                      params = params.mat,
                                                      ic = ic.mat,
                                                      data = data,
                                                      Nmc = Nmc,
                                                      draw = draw) 
pred.ic.param.process$tea <- predict_state_one_precip(type = type,
                                                    site = "Tea Control",
                                                    params = params.mat,
                                                    ic = ic.mat,
                                                    data = data,
                                                    Nmc = Nmc,
                                                    draw = draw)

conf.int.ic.param.process <- list()
conf.int.ic.param.process$green <- life.stage.ci(pred.ic.param.process$green)
conf.int.ic.param.process$henry <- life.stage.ci(pred.ic.param.process$henry)
conf.int.ic.param.process$tea <- life.stage.ci(pred.ic.param.process$tea)

all.uncertainty <- c("ic", "parameter", "process", "random effect")
pred.all <- list()
pred.all$green <- predict_state_one_precip(type = all.uncertainty,
                                         site = "Green Control",
                                         params = params.mat,
                                         ic = ic.mat,
                                         data = data,
                                         Nmc = Nmc,
                                         draw = draw) 
pred.all$henry <- predict_state_one_precip(type = all.uncertainty,
                                         site = "Henry Control",
                                         params = params.mat,
                                         ic = ic.mat,
                                         data = data,
                                         Nmc = Nmc,
                                         draw = draw) 
pred.all$tea <- predict_state_one_precip(type = all.uncertainty,
                                       site = "Tea Control",
                                       params = params.mat,
                                       ic = ic.mat,
                                       data = data,
                                       Nmc = Nmc,
                                       draw = draw) 
conf.int.all <- list()
conf.int.all$green <- life.stage.ci(pred.all$green)
conf.int.all$henry <- life.stage.ci(pred.all$henry)
conf.int.all$tea <- life.stage.ci(pred.all$tea)


# Green Control site index for data
s <- 1 
time <- 1:(data$N_est[s]-1)
time.data <- time + 1
obs.larva <- data$y[1,time.data,s]
obs.nymph <- data$y[2,time.data,s]
obs.adult <- data$y[3,time.data,s]

gg.data.larv.green <- gg_data_ci_partition(CI.param = conf.int.param$green,
                                           CI.all = conf.int.all$green,
                                           CI.ic.param.proc = conf.int.ic.param.process$green,
                                           CI.ic.param = conf.int.ic.param$green,
                                           life.stage = "Larvae",
                                           time = time,
                                           obs = obs.larva,
                                           median = conf.int.param$green$Larvae[2,])
gg.data.nymph.green <- gg_data_ci_partition(CI.param = conf.int.param$green,
                                            CI.all = conf.int.all$green,
                                            CI.ic.param.proc = conf.int.ic.param.process$green,
                                            CI.ic.param = conf.int.ic.param$green,
                                            life.stage = "Nymph",
                                            time = time,
                                            obs = obs.nymph,
                                            median = conf.int.param$green$Nymph[2,])
gg.data.adult.green <- gg_data_ci_partition(CI.param = conf.int.param$green,
                                            CI.all = conf.int.all$green,
                                            CI.ic.param.proc = conf.int.ic.param.process$green,
                                            CI.ic.param = conf.int.ic.param$green,
                                            life.stage = "Adult",
                                            time = time,
                                            obs = obs.adult,
                                            median = conf.int.param$green$Adult[2,])

r.square.green.nymph <- r_square_one2one(conf.int.ic.param.process$green$Nymph[2,], obs.nymph)
r.square.green.adult <- r_square_one2one(conf.int.ic.param.process$green$Adult[2,], obs.adult)
r.square.green.larva <- r_square_one2one(conf.int.ic.param.process$green$Larvae[2,], obs.larva)

rmse.green.nymph <- rmse(pred = conf.int.ic.param.process$green$Nymph[2,], obs = obs.nymph)
rmse.green.adult <- rmse(pred = conf.int.ic.param.process$green$Adult[2,], obs = obs.adult)
rmse.green.larvae <- rmse(pred = conf.int.ic.param.process$green$Larvae[2,], obs = obs.larva)





lm.nymph <- lm(conf.int.all$green$Nymph[2,] ~ obs.nymph)
lm.adult <- lm(conf.int.all$green$Adult[2,] ~ obs.adult)
lm.larva <- lm(conf.int.all$green$Larvae[2,] ~ obs.larva)


# Henry Control site index for data
s <- 2 
time <- 1:(data$N_est[s]-1)
time.data <- time + 1
obs.larva <- data$y[1,time.data,s]
obs.nymph <- data$y[2,time.data,s]
obs.adult <- data$y[3,time.data,s]

gg.data.larv.henry <- gg_data_ci_partition(CI.param = conf.int.param$henry,
                                           CI.all = conf.int.all$henry,
                                           CI.ic.param.proc = conf.int.ic.param.process$henry,
                                           CI.ic.param = conf.int.ic.param$henry,
                                           life.stage = "Larvae",
                                           time = time,
                                           obs = obs.larva,
                                           median = conf.int.param$henry$Larvae[2,])
gg.data.nymph.henry <- gg_data_ci_partition(CI.param = conf.int.param$henry,
                                            CI.all = conf.int.all$henry,
                                            CI.ic.param.proc = conf.int.ic.param.process$henry,
                                            CI.ic.param = conf.int.ic.param$henry,
                                            life.stage = "Nymph",
                                            time = time,
                                            obs = obs.nymph,
                                            median = conf.int.ic.param.process$henry$Nymph[2,])
gg.data.adult.henry <- gg_data_ci_partition(CI.param = conf.int.param$henry,
                                            CI.all = conf.int.all$henry,
                                            CI.ic.param.proc = conf.int.ic.param.process$henry,
                                            CI.ic.param = conf.int.ic.param$henry,
                                            life.stage = "Adult",
                                            time = time,
                                            obs = obs.adult,
                                            median = conf.int.param$henry$Adult[2,])

r.square.henry.nymph <- r_square_one2one(conf.int.ic.param.process$henry$Nymph[2,], obs.nymph)
r.square.henry.adult <- r_square_one2one(conf.int.ic.param.process$henry$Adult[2,], obs.adult)
r.square.henry.larva <- r_square_one2one(conf.int.ic.param.process$henry$Larvae[2,], obs.larva)

rmse.henry.nymph <- rmse(pred = conf.int.ic.param.process$henry$Nymph[2,], obs = obs.nymph)
rmse.henry.adult <- rmse(pred = conf.int.ic.param.process$henry$Adult[2,], obs = obs.adult)
rmse.henry.larvae <- rmse(pred = conf.int.ic.param.process$henry$Larvae[2,], obs = obs.larva)


# Tea Control site index for data
s <- 3 
time <- 1:(data$N_est[s]-1)
time.data <- time + 1
obs.larva <- data$y[1,time.data,s]
obs.nymph <- data$y[2,time.data,s]
obs.adult <- data$y[3,time.data,s]

gg.data.larv.tea <- gg_data_ci_partition(CI.param = conf.int.param$tea,
                                         CI.all = conf.int.all$tea,
                                         CI.ic.param.proc = conf.int.ic.param.process$tea,
                                         CI.ic.param = conf.int.ic.param$tea,
                                         life.stage = "Larvae",
                                         time = time,
                                         obs = obs.larva,
                                         median = conf.int.param$tea$Larvae[2,])
gg.data.nymph.tea <- gg_data_ci_partition(CI.param = conf.int.param$tea,
                                          CI.all = conf.int.all$tea,
                                          CI.ic.param.proc = conf.int.ic.param.process$tea,
                                          CI.ic.param = conf.int.ic.param$tea,
                                          life.stage = "Nymph",
                                          time = time,
                                          obs = obs.nymph,
                                          median = conf.int.param$tea$Nymph[2,])
gg.data.adult.tea <- gg_data_ci_partition(CI.param = conf.int.param$tea,
                                          CI.all = conf.int.all$tea,
                                          CI.ic.param.proc = conf.int.ic.param.process$tea,
                                          CI.ic.param = conf.int.ic.param$tea,
                                          life.stage = "Adult",
                                          time = time,
                                          obs = obs.adult,
                                          median = conf.int.param$tea$Adult[2,])

r.square.tea.nymph <- r_square_one2one(conf.int.ic.param.process$tea$Nymph[2,], obs.nymph)
r.square.tea.adult <- r_square_one2one(conf.int.ic.param.process$tea$Adult[2,], obs.adult)
r.square.tea.larva <- r_square_one2one(conf.int.ic.param.process$tea$Larvae[2,], obs.larva)

rmse.tea.nymph <- rmse(pred = conf.int.ic.param.process$tea$Nymph[2,], obs = obs.nymph)
rmse.tea.adult <- rmse(pred = conf.int.ic.param.process$tea$Adult[2,], obs = obs.adult)
rmse.tea.larvae <- rmse(pred = conf.int.ic.param.process$tea$Larvae[2,], obs = obs.larva)

plot_uncertainty <- function(data, life.stage){
  ggplot <- ggplot(data, aes(time, median)) +
    
    # null CI
    geom_ribbon(aes(ymin = ci.all.low, 
                    ymax = ci.all.high,
                    fill = "Random Effect + Process + Parameter + IC"),
                alpha = 0.8) +
    
    geom_ribbon(aes(ymin = ci.ic.param.proc.low, 
                    ymax = ci.ic.param.proc.high,
                    fill = "Process + Parameter + IC"),
                alpha = 0.8) +
    
    geom_ribbon(aes(ymin = ci.ic.param.low, 
                    ymax = ci.ic.param.high,
                    fill = "Parameter + IC"),
                alpha = 0.8) +
    
    geom_ribbon(aes(ymin = ci.param.low,
                    ymax = ci.param.high,
                    fill = "Parameter"),
                alpha = 0.8) +
    
    # median temp line
    geom_line() +
    
    # add observations
    geom_point(aes(time, obs)) +
    
    # fill geom_ribbons, also maps to legend
    scale_fill_manual(name = "Uncertainty", 
                      values = c("Random Effect + Process + Parameter + IC" = "lightblue",
                                 "Parameter + IC" = "lightgreen",
                                 "Parameter" = "purple",
                                 "Process + Parameter + IC" = "grey")) +
    
    # labels
    labs(x = "Time",
         y = "Individuals",
         title = life.stage) +
    
    # themes
    theme_classic() +
    theme(legend.position = c(0.85, 0.75), # put legent inside plot
          text = element_text(size = 20))  # increase font size 
  return(ggplot)
}


# Green Control


# green plots

plot_uncertainty(gg.data.larv.green, "Larvae") + coord_cartesian(ylim = c(0, 2000))
plot_uncertainty(gg.data.nymph.green, "Nymph") + coord_cartesian(ylim = c(0, 200))

#gg.data.adult.green$median <- conf.int.ic.param.process$green$Adult[2,]
plot_uncertainty(gg.data.adult.green, "Adult") + coord_cartesian(ylim = c(0, 25))


# Henry Control


# henry plots

plot_uncertainty(gg.data.larv.henry, "Larvae") + coord_cartesian(ylim = c(0, 2000))
plot_uncertainty(gg.data.nymph.henry, "Nymph") + coord_cartesian(ylim = c(0, 200))
plot_uncertainty(gg.data.adult.henry, "Adult") + coord_cartesian(ylim = c(0, 25)) 


# Tea Control


# Tea plots

plot_uncertainty(gg.data.larv.tea, "Larvae") + coord_cartesian(ylim = c(0, 2000))
plot_uncertainty(gg.data.nymph.tea, "Nymph") + coord_cartesian(ylim = c(0, 400))
plot_uncertainty(gg.data.adult.tea, "Adult") + coord_cartesian(ylim = c(0, 100))


## variance calculations
var.calc <- function(pred){
  var.param <- list()
  var.param$green <- apply(pred.param$green, 2, var)
  var.param$henry <- apply(pred.param$henry, 2, var)
  var.param$tea <- apply(pred.param$tea, 2, var)
  return(var.param)
}
var.param <- var.calc(pred.param)
var.param.ic <- var.calc(pred.ic.param)
var.param.ic.proc <- var.calc(pred.ic.param.process)
var.all <- var.calc(pred.all)

var.mat <- list(var.green = rbind(var.param$green, var.param.ic$green, var.param.ic.proc$green, var.all$green),
                var.henry = rbind(var.param$henry, var.param.ic$henry, var.param.ic.proc$henry, var.all$henry),
                var.tea = rbind(var.param$tea, var.param.ic$tea, var.param.ic.proc$tea, var.all$tea))









