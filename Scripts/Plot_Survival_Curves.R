library(ecoforecastR)
library(ggplot2)

load("/projectnb/dietzelab/fosterj/FinalOut/SurvivalModels/Flat_Larvae_MaxTemp_summer.RData")
flat.larv.temp <- jags.out
summary.flat.larv.temp <- summary(flat.larv.temp)

load("/projectnb/dietzelab/fosterj/FinalOut/SurvivalModels/Flat_Larvae_NULL.RData")
flat.larv.null <- jags.out
summary.flat.larv.null <- summary(flat.larv.null)

load("/projectnb/dietzelab/fosterj/FinalOut/SurvivalModels/Flat_Larvae_Precip_summer.RData")
flat.larv.precip <- jags.out
summary.larv.precip <- summary(flat.larv.precip)

mat.flat.nymph.temp <- as.matrix(flat.larv.temp)
mat.flat.larv.null <- as.matrix(flat.larv.null)
mat.flat.nymph.precip <- as.matrix(flat.larv.precip)


rm(jags.out)

# survival process model
"
## loop for daily survival for each row
for(i in 1:N){
  for(t in 1:(N_days[i]-1)){
    logit(lambda[i,t]) <- lambda.mu + 
      alpha[site.index[i]] + 
      beta*met[3,i,start.index[i]+t-1]
  }
}
"

ci.null <- apply(mat.flat.larv.null, 2, quantile, c(0.025, 0.5, 0.975))

### temp curves

# theoretical temp values (deg C)
temp <- seq(-10, 30, by = 0.1)

Nmc <- 1000
draw <- sample(nrow(mat.flat.nymph.temp), Nmc, replace = TRUE)

surv <- array(data = NA, dim = c(Nmc, length(temp), 2))
theta <- mat.flat.nymph.temp[draw,]
for(s in 2:3){
  alpha <- theta[,paste("alpha[", s, "]", sep = "")]
  for(t in 1:length(temp)){
    surv[,t,s-1] <- inv.logit(theta[,"lambda.mu"] +
                      theta[,"beta"]*temp[t] +
                      alpha)
  }
}

ci.site.1 <- apply(surv[,,1], 2, quantile, c(0.025, 0.5, 0.975))
ci.site.2 <- apply(surv[,,2], 2, quantile, c(0.025, 0.5, 0.975))

plot(temp, ci.site.1[2,])

gg.data <- data.frame(null.low = rep(inv.logit(ci.null[1,"lambda.mu"]), ncol(ci.site.1)),
                      null.med = rep(inv.logit(ci.null[2,"lambda.mu"]), ncol(ci.site.1)),
                      null.up = rep(inv.logit(ci.null[3,"lambda.mu"]), ncol(ci.site.1)),
                      ci.temp.low = ci.site.1[1,],
                      ci.temp.med = ci.site.1[2,],
                      ci.temp.up = ci.site.1[3,],
                      temp = temp)
alpha = 0.8
col.null <- "#67a9cf"
col.met <- "#ef8a62"
ggplot <- ggplot(gg.data, aes(temp, ci.temp.med)) +
  
  # null CI
  geom_ribbon(aes(ymin = null.low, 
                  ymax = null.up,
                  fill = "Null"),
              alpha = alpha) +
  
  # temp curve CI
  geom_ribbon(aes(ymin = ci.temp.low, 
                  ymax = ci.temp.up,
                  fill = "Temperature"),
              alpha = alpha) +
  
  # median temp line
  geom_line() +
  
  # median null line
  geom_line(aes(temp, null.med)) +
  
  # fill geom_ribbons, also maps to legend
  scale_fill_manual(name = "", 
                    values = c(col.null, col.met)) +

  # labels
  labs(x = "Temperature (Deg C)",
       y = "Daily Survival Probablity",
       title = "Flat Larvae Survival Curve") +
  
  # themes
  theme_classic() +
  theme(legend.position = c(0.75, 0.25), # put legent inside plot
        text = element_text(size = 20))  # increase font size 
ggplot

## Precipitation Curves

temp <- seq(0, 80, by = 0.1)

Nmc <- 5000
draw <- sample(nrow(mat.flat.nymph.precip), Nmc, replace = TRUE)

surv <- array(data = NA, dim = c(Nmc, length(temp), 2))
theta <- mat.flat.nymph.precip[draw,]
for(s in 2:3){
  alpha <- theta[,paste("alpha[", s, "]", sep = "")]
  for(t in 1:length(temp)){
    surv[,t,s-1] <- inv.logit(theta[,"lambda.mu"] +
                                theta[,"beta"]*temp[t] +
                                alpha)
  }
}

ci.site.1 <- apply(surv[,,1], 2, quantile, c(0.025, 0.5, 0.975))
ci.site.2 <- apply(surv[,,2], 2, quantile, c(0.025, 0.5, 0.975))

plot(temp, ci.site.1[2,])

gg.data <- data.frame(null.low = rep(inv.logit(ci.null[1,"lambda.mu"]), ncol(ci.site.1)),
                      null.med = rep(inv.logit(ci.null[2,"lambda.mu"]), ncol(ci.site.1)),
                      null.up = rep(inv.logit(ci.null[3,"lambda.mu"]), ncol(ci.site.1)),
                      ci.temp.low = ci.site.1[1,],
                      ci.temp.med = ci.site.1[2,],
                      ci.temp.up = ci.site.1[3,],
                      temp = temp)
alpha = 0.8
col.null <- "#67a9cf"
col.met <- "#ef8a62"
ggplot <- ggplot(gg.data, aes(temp, ci.temp.med)) +
  
  # null CI
  geom_ribbon(aes(ymin = null.low, 
                  ymax = null.up,
                  fill = "Null"),
              alpha = alpha) +
  
  # temp curve CI
  geom_ribbon(aes(ymin = ci.temp.low, 
                  ymax = ci.temp.up,
                  fill = "Precipitation"),
              alpha = alpha) +
  
  # median temp line
  geom_line() +
  
  # median null line
  geom_line(aes(temp, null.med)) +
  
  # fill geom_ribbons, also maps to legend
  scale_fill_manual(name = "", 
                    values = c(col.null, col.met)) +
  
  # labels
  labs(x = "Precipitation (mm)",
       y = "Daily Survival Probablity",
       title = "Flat Larvae Survival Curve") +
  
  # themes
  theme_classic() +
  theme(legend.position = c(0.75, 0.25), # put legent inside plot
        text = element_text(size = 20))  # increase font size 
ggplot

###########################################################################
###########                     NYMPHS                      ###############
###########################################################################

load("/projectnb/dietzelab/fosterj/FinalOut/SurvivalModels/Flat_Nymph_MaxTemp.RData")
flat.nymph.temp <- jags.out
summary.flat.nymph.temp <- summary(flat.nymph.temp)

load("/projectnb/dietzelab/fosterj/FinalOut/SurvivalModels/Flat_Nymph_NULL.RData")
flat.nymph.null <- jags.out
summary.flat.nymph.null <- summary(flat.nymph.null)

load("/projectnb/dietzelab/fosterj/FinalOut/SurvivalModels/Flat_Nymph_Precip.RData")
flat.nymph.precip <- jags.out
summary.nymph.precip <- summary(flat.nymph.precip)

mat.flat.nymph.temp <- as.matrix(flat.nymph.temp)
mat.flat.nymph.null <- as.matrix(flat.nymph.null)
mat.flat.nymph.precip <- as.matrix(flat.nymph.precip)

ci.null <- apply(mat.flat.larv.null, 2, quantile, c(0.025, 0.5, 0.975))

### temp curves

# theoretical temp values (deg C)
temp <- seq(-10, 30, by = 0.1)

Nmc <- 5000
draw <- sample(nrow(mat.flat.nymph.temp), Nmc, replace = TRUE)

surv <- array(data = NA, dim = c(Nmc, length(temp), 3))
theta <- mat.flat.nymph.temp[draw,]
for(s in 1:3){
  alpha <- theta[,paste("alpha[", s, "]", sep = "")]
  for(t in 1:length(temp)){
    surv[,t,s] <- inv.logit(theta[,"lambda.mu"] +
                                theta[,"beta"]*temp[t] +
                                alpha)
  }
}

ci.site.1 <- apply(surv[,,1], 2, quantile, c(0.025, 0.5, 0.975))
ci.site.2 <- apply(surv[,,2], 2, quantile, c(0.025, 0.5, 0.975))
ci.site.3 <- apply(surv[,,3], 2, quantile, c(0.025, 0.5, 0.975))

surv.null <- matrix(NA, Nmc, 3)
theta <- mat.flat.nymph.null[draw,]
for(s in 1:3){
  alpha <- theta[,paste("alpha[", s, "]", sep = "")]
  surv.null[,s] <- inv.logit(theta[,"lambda.mu"] + alpha)
}

ci.site.null.1 <- quantile(surv.null[,1], c(0.025, 0.5, 0.975))
ci.site.null.2 <- quantile(surv.null[,2], c(0.025, 0.5, 0.975))
ci.site.null.3 <- quantile(surv.null[,3], c(0.025, 0.5, 0.975))

gg.data <- data.frame(null.low = rep(ci.site.null.1[1], ncol(ci.site.1)),
                      null.med = rep(ci.site.null.1[2], ncol(ci.site.1)),
                      null.up = rep(ci.site.null.1[3], ncol(ci.site.1)),
                      null.mean = rep(mean(surv.null[,1]), ncol(ci.site.1)),
                      ci.temp.low = ci.site.1[1,],
                      ci.temp.med = ci.site.1[2,],
                      ci.temp.up = ci.site.1[3,],
                      temp = temp)
alpha = 0.75
col.null <- "#67a9cf"
col.met <- "#ef8a62"
ggplot <- ggplot(gg.data, aes(temp, ci.temp.med)) +
  
  # temp curve CI
  geom_ribbon(aes(ymin = ci.temp.low, 
                  ymax = ci.temp.up,
                  fill = "Temperature"),
              alpha = alpha) +
  
  # null CI
  geom_ribbon(aes(ymin = null.low, 
                  ymax = null.up,
                  fill = "Null"),
              alpha = alpha) +
  
  
  # median temp line
  geom_line() +
  
  # median null line
  geom_line(aes(temp, null.mean),
            linetype = "dashed") +
  
  # fill geom_ribbons, also maps to legend
  scale_fill_manual(name = "", 
                    values = c(col.null, col.met)) +
  
  # labels
  labs(x = "Temperature (Deg C)",
       y = "Daily Survival Rate",
       title = "Flat Nymph Survival Curve") +
  
  # themes
  theme_classic() +
  coord_cartesian(ylim = c(0.9,1)) +
  theme(legend.position = c(0.35, 0.25), # put legent inside plot
        text = element_text(size = 24))  # increase font size 
ggplot

## Precipitation Curves

temp <- seq(0, 5, by = 0.01)

draw <- sample(nrow(mat.flat.nymph.precip), Nmc, replace = TRUE)

surv <- array(data = NA, dim = c(Nmc, length(temp), 3))
theta <- mat.flat.nymph.precip[draw,]
for(s in 1:3){
  alpha <- theta[,paste("alpha[", s, "]", sep = "")]
  for(t in 1:length(temp)){
    surv[,t,s] <- inv.logit(theta[,"lambda.mu"] +
                                theta[,"beta"]*temp[t] +
                                alpha)
  }
}

ci.site.1 <- apply(surv[,,1], 2, quantile, c(0.025, 0.5, 0.975))
ci.site.2 <- apply(surv[,,2], 2, quantile, c(0.025, 0.5, 0.975))

surv.null <- matrix(NA, Nmc, 3)
theta <- mat.flat.nymph.null[draw,]
for(s in 1:3){
  alpha <- theta[,paste("alpha[", s, "]", sep = "")]
  surv.null[,s] <- inv.logit(theta[,"lambda.mu"] + alpha)
}

ci.site.null.1 <- quantile(surv.null[,1], c(0.025, 0.5, 0.975))
ci.site.null.2 <- quantile(surv.null[,2], c(0.025, 0.5, 0.975))
ci.site.null.3 <- quantile(surv.null[,3], c(0.025, 0.5, 0.975))

gg.data <- data.frame(null.low = rep(ci.site.null.1[1], ncol(ci.site.1)),
                      null.med = rep(ci.site.null.1[2], ncol(ci.site.1)),
                      null.up = rep(ci.site.null.1[3], ncol(ci.site.1)),
                      null.mean = rep(mean(surv.null[,1]), ncol(ci.site.1)),
                      ci.temp.low = ci.site.1[1,],
                      ci.temp.med = ci.site.1[2,],
                      ci.temp.up = ci.site.1[3,],
                      temp = temp)

alpha <- 0.75
col.null <- "#67a9cf"
col.met <- "#ef8a62"
ggplot <- ggplot(gg.data, aes(temp, ci.temp.med)) +
  
  # null CI
  geom_ribbon(aes(ymin = null.low, 
                  ymax = null.up,
                  fill = "Null"),
              alpha = alpha) +
  
  # temp curve CI
  geom_ribbon(aes(ymin = ci.temp.low, 
                  ymax = ci.temp.up,
                  fill = "Precipitation"),
              alpha = alpha) +
  
  # median temp line
  geom_line() +
  
  # median null line
  geom_line(aes(temp, null.mean),
            linetype = "dashed") +
  
  # fill geom_ribbons, also maps to legend
  scale_fill_manual(name = "", 
                    values = c(col.null, col.met)) +
  
  # labels
  labs(x = "Precipitation (mm)",
       y = "Daily Survival Rate",
       title = "Flat Nymph Survival Curve") +
  
  # themes
  theme_classic() +
  coord_cartesian(ylim = c(0.9,1)) +
  theme(legend.position = c(0.65, 0.25), # put legent inside plot
        text = element_text(size = 24))  # increase font size 
ggplot













