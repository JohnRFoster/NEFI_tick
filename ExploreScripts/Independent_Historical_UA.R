
library(knitr)
library(ecoforecastR)
library(runjags)
library(fosteR)
library(mvtnorm)
library(ggplot2)
library(grid)
library(gridExtra)
library(boot)

source("Functions/Future_Met.R")
source("Functions/life_stage_ci.R")
source("Functions/model_performance_metrics.R")
source("Functions/ammend_chains.R")
source("Functions/cary_tick_met_JAGS.R")
source("Functions/plot_post_prior.R")
source("Functions/site_data_met.R")
source("Functions/gg_data_ci_partition.R")
source("Models/predict_one_null.R")

dir <- "../FinalOut/A_Correct/NULL/"
met.variable <- NULL

params.ls <- predict.ls <- list()
sites <- c("Green Control", "Henry Control", "Tea Control")
pattern <- " Control"
for(i in 1:3){
  model <- "Combined_thinMat_NULL_"
  site <- sites[i]
  site.folder <- gsub(pattern, "", site)
  model <- paste0(model, site.folder, "Control.RData")
  load(file = file.path(dir, site.folder, model))
  params.ls[[i]] <- params.mat
  predict.ls[[i]] <- predict.mat
}


Nmc <- 300
draw <- sample.int(nrow(params.ls[[1]]), Nmc, replace = TRUE)
data <- cary_ticks_met_JAGS()


pred.param <- list()
conf.int.param <- list()
for(i in 1:3){
  pred.param[[i]] <- predict_one_null(type = "parameter",
                                                         site = sites[i],
                                                         # met.variable = met.variable,
                                                         params = params.ls[[i]],
                                                         ic = predict.ls[[i]],
                                                         data = data,
                                                         Nmc = Nmc,
                                                         draw = draw)
  conf.int.param[[i]] <- life.stage.ci(pred.param[[i]][[1]])
}

pred.ic <- list()
conf.int.ic <- list()
for(i in 1:3){
  pred.ic[[i]] <- predict_one_null(type = "ic",
                                                         site = sites[i],
                                                         # met.variable = met.variable,
                                                         params = params.ls[[i]],
                                                         ic = predict.ls[[i]],
                                                         data = data,
                                                         Nmc = Nmc,
                                                         draw = draw)
  conf.int.ic[[i]] <- life.stage.ci(pred.ic[[i]][[1]])
}


# param + IC
type <- c("parameter", "ic")
pred.ic.param <- list()
conf.int.ic.param <- list()
for(i in 1:3){
  pred.ic.param[[i]] <- predict_one_null(type = type,
                                                     site = sites[i],
                                                     # met.variable = met.variable,
                                                     params = params.ls[[i]],
                                                     ic = predict.ls[[i]],
                                                     data = data,
                                                     Nmc = Nmc,
                                                     draw = draw)
  conf.int.ic.param[[i]] <- life.stage.ci(pred.ic.param[[i]][[1]])
}

# param + IC + Process
type <- c("parameter", "ic", "process")
pred.ic.param.process <- list()
conf.int.ic.param.process <- list()
for(i in 1:3){
  pred.ic.param.process[[i]] <- predict_one_null(type = type,
                                                           site = sites[i],
                                                           # met.variable = met.variable,
                                                           params = params.ls[[i]],
                                                           ic = predict.ls[[i]],
                                                           data = data,
                                                           Nmc = Nmc,
                                                           draw = draw)
  conf.int.ic.param.process[[i]] <- life.stage.ci(pred.ic.param.process[[i]][[1]])
}

# all uncertainty
# type <- c("ic", "parameter", "process", "random effect")
# pred.ic.param.process.re <- list()
# conf.int.ic.param.process.re <- list()
# for(i in 1:3){
#   pred.ic.param.process.re[[i]] <- predict_one_null(type = type,
#                                                                    site = sites[i],
#                                                                    met.variable = met.variable,
#                                                                    params = params.ls[[i]],
#                                                                    ic = predict.ls[[i]],
#                                                                    data = data,
#                                                                    Nmc = Nmc,
#                                                                    draw = draw)
#   conf.int.ic.param.process.re[[i]] <- life.stage.ci(pred.ic.param.process.re[[i]][[1]])
# }
ls.names <- c("green", "henry", "tea")
names(conf.int.ic.param) <- ls.names
names(conf.int.ic) <- ls.names
names(conf.int.param) <- ls.names
names(conf.int.ic.param.process) <- ls.names

# Green Control site index for data
s <- 1 
time <- 1:(data$N_est[s]-1)
time.data <- time + 1
obs.larva <- data$y[1,time.data,s]
obs.nymph <- data$y[2,time.data,s]
obs.adult <- data$y[3,time.data,s]

gg.data.larv.green <- gg_data_ci_partition(CI.param = conf.int.param$green,
                                           CI.all = NULL,
                                           CI.ic.param.proc = conf.int.ic.param.process$green,
                                           CI.ic.param = conf.int.ic.param$green,
                                           life.stage = "Larvae",
                                           time = time,
                                           obs = obs.larva,
                                           median = conf.int.param$green$Larvae[2,])
gg.data.nymph.green <- gg_data_ci_partition(CI.param = conf.int.param$green,
                                            CI.all = NULL,
                                            CI.ic.param.proc = conf.int.ic.param.process$green,
                                            CI.ic.param = conf.int.ic.param$green,
                                            life.stage = "Nymph",
                                            time = time,
                                            obs = obs.nymph,
                                            median = conf.int.param$green$Nymph[2,])
gg.data.adult.green <- gg_data_ci_partition(CI.param = conf.int.param$green,
                                            CI.all = NULL,
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
                                           CI.all = NULL,
                                           CI.ic.param.proc = conf.int.ic.param.process$henry,
                                           CI.ic.param = conf.int.ic.param$henry,
                                           life.stage = "Larvae",
                                           time = time,
                                           obs = obs.larva,
                                           median = conf.int.param$henry$Larvae[2,])
gg.data.nymph.henry <- gg_data_ci_partition(CI.param = conf.int.param$henry,
                                            CI.all = NULL,
                                            CI.ic.param.proc = conf.int.ic.param.process$henry,
                                            CI.ic.param = conf.int.ic.param$henry,
                                            life.stage = "Nymph",
                                            time = time,
                                            obs = obs.nymph,
                                            median = conf.int.ic.param.process$henry$Nymph[2,])
gg.data.adult.henry <- gg_data_ci_partition(CI.param = conf.int.param$henry,
                                            CI.all = NULL,
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
                                         CI.all = NULL,
                                         CI.ic.param.proc = conf.int.ic.param.process$tea,
                                         CI.ic.param = conf.int.ic.param$tea,
                                         life.stage = "Larvae",
                                         time = time,
                                         obs = obs.larva,
                                         median = conf.int.param$tea$Larvae[2,])
gg.data.nymph.tea <- gg_data_ci_partition(CI.param = conf.int.param$tea,
                                          CI.all = NULL,
                                          CI.ic.param.proc = conf.int.ic.param.process$tea,
                                          CI.ic.param = conf.int.ic.param$tea,
                                          life.stage = "Nymph",
                                          time = time,
                                          obs = obs.nymph,
                                          median = conf.int.param$tea$Nymph[2,])
gg.data.adult.tea <- gg_data_ci_partition(CI.param = conf.int.param$tea,
                                          CI.all = NULL,
                                          CI.ic.param.proc = conf.int.ic.param.process$tea,
                                          CI.ic.param = conf.int.ic.param$tea,
                                          life.stage = "Adult",
                                          time = time,
                                          obs = obs.adult,
                                          median = conf.int.param$tea$Adult[2,])

r.square.tea.nymph <- r_square_one2one(conf.int.ic.param.process$tea$Nymph[2,], obs.nymph)
r.square.tea.adult <- r_square_one2one(conf.int.ic.param.process$tea$Adult[2,], obs.adult)
r.square.tea.larva <- r_square_one2one(conf.int.ic.param.process$tea$Larvae[2,], obs.larva)



rmse <- function(pred, obs){
  n <- length(pred)
  rmse <- sqrt((1/n)*((pred - obs)^2))
  return(rmse)
}


rmse.tea.nymph <- rmse(pred = conf.int.ic.param.process$tea$Nymph[2,], obs = obs.nymph)
rmse.tea.adult <- rmse(pred = conf.int.ic.param.process$tea$Adult[2,], obs = obs.adult)
rmse.tea.larvae <- rmse(pred = conf.int.ic.param.process$tea$Larvae[2,], obs = obs.larva)

plot_uncertainty <- function(data, life.stage){
  ggplot <- ggplot(data, aes(time, median)) +
    
    # null CI
    # geom_ribbon(aes(ymin = ci.all.low, 
    #                 ymax = ci.all.high,
    #                 fill = "Random Effect + Process + Parameter + IC"),
    #             alpha = 0.8) +
    
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
    theme(legend.position = c(0.25, 0.85), # put legent inside plot
          text = element_text(size = 20))  # increase font size 
  return(ggplot)
}


# Green Control


# green plots

plot_uncertainty(gg.data.larv.green, "Larvae") + coord_cartesian(ylim = c(0, 2000))
plot_uncertainty(gg.data.nymph.green, "Nymph") + coord_cartesian(ylim = c(0, 100))

#gg.data.adult.green$median <- conf.int.ic.param.process$green$Adult[2,]
plot_uncertainty(gg.data.adult.green, "Adult") + coord_cartesian(ylim = c(0, 25))


# Henry Control


# henry plots

plot_uncertainty(gg.data.larv.henry, "Larvae") + coord_cartesian(ylim = c(0, 2000))
plot_uncertainty(gg.data.nymph.henry, "Nymph") + coord_cartesian(ylim = c(0, 100))
plot_uncertainty(gg.data.adult.henry, "Adult") + coord_cartesian(ylim = c(0, 25)) 


# Tea Control


# Tea plots

plot_uncertainty(gg.data.larv.tea, "Larvae") + coord_cartesian(ylim = c(0, 2000))
plot_uncertainty(gg.data.nymph.tea, "Nymph") + coord_cartesian(ylim = c(0, 200))
plot_uncertainty(gg.data.adult.tea, "Adult") + coord_cartesian(ylim = c(0, 100))


## variance calculations
var.calc <- function(pred){
  var.param <- list()
  var.param$green <- apply(pred[1,,], 1, var)
  var.param$henry <- apply(pred, 2, var)
  var.param$tea <- apply(pred, 2, var)
  return(var.param)
}

var.param <- var.calc(pred.param)
var.param.ic <- var.calc(pred.ic.param)
var.param.ic.proc <- var.calc(pred.ic.param.process)
var.all <- var.calc(pred.all)

var.mat <- list(var.green = rbind(var.param$green, var.param.ic$green, var.param.ic.proc$green, var.all$green),
                var.henry = rbind(var.param$henry, var.param.ic$henry, var.param.ic.proc$henry, var.all$henry),
                var.tea = rbind(var.param$tea, var.param.ic$tea, var.param.ic.proc$tea, var.all$tea))

var.param <- var.param.ic <- var.param.ic.proc <- list()
for(ls in 1:3){
  var.param[[ls]] <- apply(pred.ic[[2]][ls,,], 1, var)
  var.param.ic[[ls]] <- apply(pred.ic.param[[2]][ls,,], 1, var)
  var.param.ic.proc[[ls]] <- apply(pred.ic.param.process[[2]][ls,,], 1, var)
}
var.mat <- list()
for(ls in 1:3){
  var.mat[[ls]] <- rbind(var.param[[ls]], var.param.ic[[ls]], var.param.ic.proc[[ls]])  
}
time.1 <- 1:72
V.pred.rel <- apply(var.mat[[2]],2,function(x) {x/max(x)})
# V.pred.rel <- apply(V.pred.rel, 2, sort)
plot(time.1,V.pred.rel[1,],ylim=c(0,1),type='n',main="Relative Variance: In-Sample (Nymph)",ylab="Proportion of Variance",xlab="time")
ciEnvelope(time.1, rep(0,ncol(V.pred.rel)),V.pred.rel[1,],col="purple")
ciEnvelope(time.1, V.pred.rel[1,],V.pred.rel[2,],col="lightgreen")
ciEnvelope(time.1, V.pred.rel[2,],V.pred.rel[3,],col="grey")


plot(conf.int.ic.param.process$tea$Nymph[2,], obs.nymph)
plot(conf.int.ic.param.process$tea$Adult[2,], obs.adult)
plot(conf.int.ic.param.process$tea$Nymph[2,], obs.nymph)







