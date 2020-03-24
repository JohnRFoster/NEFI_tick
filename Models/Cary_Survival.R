###-------------------------------------------------------###
### Survival HB sub model for Cary ticks                  ### 
### Estimates daily survival for each site                ###
###                                                       ###
### Subset: Flat Nymphs                                   ###
### Drivers: Max Temp                                     ###
###-------------------------------------------------------###

library(ecoforecastR)
source("Functions/surv_type_iButton.R")
source("Functions/RunMCMC_Model.R")

## file path to output folder
out.folder <- "../FinalOut/SurvivalModels/Flat_Nymph/"

## set the subset of ticks to estimate survival
subset <- "Flat Nymph"
subset.out <- gsub(pattern = " ", replacement = "_", x = subset)

## set season (only need for larvae)
# season <- "winter"

# set driver
# driver.vec <- c("temp", "rh", "vpd")
# xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# driver <- driver.vec[xx]
position <- "Above"
# out.name <- paste(subset.out,driver,position, "AllTicks_NoRandom",sep = "_")

# null run
driver <- "NULL"
out.name <- paste(subset.out,driver, "AllTicks_NoRandom",sep = "_")

out.path <- paste(out.folder,out.name, sep = "")

# read data
data <- surv_type_iButton(subset, position)
# data$b_mu <- 0
# data$b_prec <- 0.01
data$b_mu <- 0
data$b_prec <- 0.001

if(driver == "temp"){
  data$met <- data$temp
} else if(driver == "rh"){
  data$met <- data$rh
} else {
  data$met <- data$vpd
} 
data[c("temp", "rh", "vpd")] <- NULL

cat("Number of cores:", nrow(data$y), "\n")

# force summary stats to be calculated
# runjags.options(force.summary = TRUE)

## met array dimensions:
  # dim 1 = met variable; 
    # 1 = Tmax
    # 2 = Tmin
    # 3 = Precip
  # dim 2 = soil cores
  # dim 3 = date

survival.model <- " model {

## priors
lambda.mu ~ dnorm(0,0.001) # global hazard probability 
tau ~ dgamma(0.01, 0.01)
beta ~ dnorm(b_mu, b_prec)

## site random effect
#alpha[1] <- 0
#for(i in 1:3){
 # alpha[i] ~ dnorm(0, tau) ## random site effect
#}

## loop for daily survival for each row
for(i in 1:N){
  for(t in 2:(N_days[i])){
    logit(lambda[i,t]) <- beta
                          # beta[1] + 
                          # alpha[site.index[i]] + 
                          # beta[2]*met[t,i] +
                          # beta[3]*(met[t,i] - met[t-1,i])
  }
}


## loop for aggrigated survival across time in each row
for(t in 1:N){
  lambda[t,1] ~ dbeta(20,1)
  log(phi[t]) <- sum(log(lambda[t,1:(N_days[t]-1)]))
}

## Binomial process for ticks surviving each day
for(i in 1:N){
  y[i,2] ~ dbinom(phi[i], y[i,1])
}

}" # end model

inits <- function(){list(beta = rnorm(3,0,0.001),
                       #  alpha = rnorm(3,0,0.1),
             		         tau = rgamma(1,0.01,0.01))}
             		       #  lambda.mu = rnorm(1,5,0.5))}

monitor <- c(#"alpha",
             #"lambda.mu",
             #"lambda",
             # "tau",
             "beta")

load.module("glm")
j.model <- jags.model(file = textConnection(survival.model),
                      data = data,
                    #  inits = inits,
                      n.adapt = 50000,
                      n.chains = 3)

jags.out <- runMCMC_Model(j.model = j.model, variableNames = monitor)

save(jags.out, j.model, file = paste(out.path, ".RData", sep = ""))

print("-------------- DIC ----------------")

cat("\nDIC for", subset, driver, "\n")
dic.samples(j.model, 10000)

print("-----------------------------------------------------------------")
print("                           END SCRIPT                            ")
print("-----------------------------------------------------------------")



#summary(run.jags)
#summary(vars.sum)
## split output
#out <- list(params = NULL, predict = NULL)
#mfit <- as.matrix(jags.out, chains = TRUE)
#pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
#chain.col <- which(colnames(mfit) == "CHAIN")
#out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
#out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])

#summary(out$params)
#summary(out$predict)
