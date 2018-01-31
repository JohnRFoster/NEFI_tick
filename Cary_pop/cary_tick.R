setwd("C:/Users/foste/Desktop/R_work/R/Cary")
library(tidyverse)
library(lubridate)
library(rjags)
library(coda)


gflu = read.csv("http://www.google.org/flutrends/about/data/flu/us/data.txt",skip=11)
time = as.Date(gflu$Date)
y = gflu$Massachusetts




data <- read.csv("Mouse_mast_Tick_Drag.csv", header = TRUE)
tick.df <- data

tick.df$Date <- as.character(tick.df$Date)
tick.df[135, 3] <- "6/12/02"                  # these three rows have two sampling dates (6/12/02-6/13/02)
tick.df[290, 3] <- "6/23/03"                  # coerced to just first day of the sampling effort 
tick.df[445, 3] <- "10/25/04"
tick.df$Date <- mdy(tick.df$Date)
tick.df <- tick.df %>% mutate(day.of.year = yday(Date)) # add day-of-year column to data frame

# data frames for each control site
green.ctrl <- tick.df %>% filter(Grid == "Green Control")    
henry.ctrl <- tick.df %>% filter(Grid == "Henry Control")
tea.ctrl <- tick.df %>% filter(Grid == "Tea Control")

# date vectors for each site
time.green <- green.ctrl[, "Date"]
time.henry <- henry.ctrl[, "Date"]
time.tea <- tea.ctrl[, "Date"]

# nymph density vectors for each control site (y in data model)
nymph.green <- green.ctrl[, "Nymphs.m2"]
nymph.henry <- henry.ctrl[, "Nymphs.m2"]
nymph.tea <- tea.ctrl[, "Nymphs.m2"]

########################################################
RandomWalk = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dnorm(x[i],tau_obs)
}

#### Process Model
for(i in 2:n){
x[i]~dnorm(x[i-1],tau_add)
}

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
}
"
### Green Control ###
data.green <- list(y = nymph.green, n = length(nymph.green), x_ic = 0.06, tau_ic = 1/0.06, a_obs = 1, 
               r_obs = 0.01, a_add = 1, r_add = 0.01)

# I chose x_ic = 0.06 as it appeared to be the median value from Allan et al. 2003 


nchain = 5
init.green <- list()
for(i in 1:nchain){
  y.samp = sample(nymph.green, length(nymph.green), replace=TRUE)
  init.green[[i]] <- list(tau_add = 1/var(y.samp) ,tau_obs = 5/var(y.samp))
}

j.model.green <- jags.model(file = textConnection(RandomWalk),
                         data = data.green,
                         inits = init.green,
                         n.chains = 5)

jags.out.green <- coda.samples(model = j.model.green,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out.green)

# Run again, but add 'x' output and increase number of iterations to 10000
jags.out.green <- coda.samples(model = j.model.green,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)

time.rng = c(1,length(time.green)) ## adjust to zoom in and out
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
out.green <- as.matrix(jags.out.green)
ci <- apply(out.green[,3:ncol(out.green)],2,quantile,c(0.025,0.5,0.975))

par(mfrow = c(1,1))
plot(time.green,ci[2,],type='n',ylim=range(nymph.green,na.rm=TRUE),ylab="Nymph Density m^-2",
     xlim=time.green[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time.green[time.rng[1]],time.green[time.rng[2]],by='month'), format = "%Y-%m")
}
ciEnvelope(time.green,ci[1,],ci[3,],col="lightBlue")
points(time.green,nymph.green,pch="+",cex=0.5)


layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
hist(1/sqrt(out.green[,1]),main=colnames(out.green)[1])
hist(1/sqrt(out.green[,2]),main=colnames(out.green)[2])
plot(out.green[,1],out.green[,2],pch=".",xlab=colnames(out.green)[1],ylab=colnames(out.green)[2])
cor(out.green[,1:2])
