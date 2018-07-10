library(rjags)
library(tidyverse)
t <- read.csv("Cary_pop/tick_cleaned")

t <- subset(t, Grid == "Green Control")

ymd <- as.data.frame(strsplit(as.character(t$DATE),"-"))
year.index <- t(ymd[1,])
year.index <- apply(year.index, 2, as.integer)

df <- as.Date(t$DATE)
df <- diff(df)
df <- c(1, df)
dt.index <- cumsum(as.integer(df))

t <- t[,c("n_larvae", "n_nymphs", "n_adults")]
dat.t <- t(t)

### average annual
year <- 1995:2005 ## years we want met
met <- read.csv("Cary_pop/Met_Cary") %>% 
  separate(DATE, c("year", "month", "day"), sep="-") %>% 
  filter(year == year) %>% 
  select(year,MAX_TEMP,MAX_RH,TOT_PREC)
met$year <- as.integer(met$year)

met <- read.csv("Cary_pop/Met_Cary")
met <- met[1097:nrow(met),c("DATE","MAX_TEMP","MAX_RH","TOT_PREC")]
met.y <- as.data.frame(strsplit(as.character(met$DATE),"-"))
met.y <- t(met.y[1,])
met.y <- apply(met.y, 2, as.integer)
met <- cbind(met.y, met[,-1])
colnames(met) <- c("year","MAX_TEMP","MAX_RH","TOT_PREC")

annual.met <- data.frame()
for(i in 1:length(year)){
  xx <- met[met$year==year[i],]
  for(c in 1:ncol(met)){
    annual.met[i, c] <- mean(xx[,c], na.rm = TRUE)
  }
}
colnames(annual.met) <- colnames(met)

precip <- annual.met[, "TOT_PREC"] ## mean annual precip
temp <- annual.met[, "MAX_TEMP"] ## mean annual temperature
rh <- annual.met[, "MAX_RH"] ## mean annual relative humidity

R <- diag(1, 3, 3)

data <- list(y = dat.t, 
             R = R,
             time = ncol(dat.t),
             year = year,
             year.index = year.index,
             temp = temp,
             precip = precip,
             rh = rh,
             b0 = as.vector(c(0,0,0)),   # regression beta means
             Vb = solve(diag(10000,3)))   # regression beta precisions

inits <- function(){list(x = dat.t)}

tick.met = " 
model {

# priors
phi.l.mu ~ dnorm(0, 10)       # larvae survival
phi.n.mu ~ dnorm(0, 10)       # nymph survival
phi.a.mu ~ dnorm(0, 10)       # adult survival
grow.ln.mu ~ dnorm(0, 10)     # larvae -> nymph transition (*truncated(phi.l))
grow.na.mu ~ dnorm(0, 10)     # nymph -> adult transition
repro.mu ~ dnorm(0, 10)       # adult -> larvae transition (reproduction)
SIGMA ~ dwish(R, 3)           # variance matrix for mvn [3 x 3] default SIGMA will be *precision*
beta.pl ~ dmnorm(b0,Vb)       # prior regression params larvae survival
beta.pn ~ dmnorm(b0,Vb)       # prior regression params nymph survival
beta.pa ~ dmnorm(b0,Vb)       # prior regression params adult survival
beta.gln ~ dmnorm(b0,Vb)      # prior regression params larvae -> nymph transition
beta.gna ~ dmnorm(b0,Vb)      # prior regression params nymph -> adult transition
beta.r ~ dmnorm(b0,Vb)        # prior regression params reproduction

for(t in 1:11) {
  logit(phi.l[t]) <- phi.l.mu + beta.pl[1]*precip[t] + beta.pl[2]*temp[t] + beta.pl[3]*rh[t] 
  logit(phi.n[t]) <- phi.n.mu + beta.pn[1]*precip[t] + beta.pn[2]*temp[t] + beta.pn[3]*rh[t] 
  logit(phi.a[t]) <- phi.a.mu + beta.pa[1]*precip[t] + beta.pa[2]*temp[t] + beta.pa[3]*rh[t] 
  logit(grow.ln[t]) <- grow.ln.mu + beta.gln[1]*precip[t] + beta.gln[2]*temp[t] + beta.gln[3]*rh[t] 
  logit(grow.na[t]) <- grow.na.mu + beta.gna[1]*precip[t] + beta.gna[2]*temp[t] + beta.gna[3]*rh[t]
  log(repro[t]) <- repro.mu + beta.r[1]*precip[t] + beta.r[2]*temp[t] + beta.r[3]*rh[t] 
}

# define transition matrix 
for(t in 1:11){
  A[1, 1, year[t]] <- phi.l[t]
  A[1, 2, year[t]] <- 0
  A[1, 3, year[t]] <- repro[t]
  A[2, 1, year[t]] <- grow.ln[t]
  A[2, 2, year[t]] <- phi.n[t]
  A[2, 3, year[t]] <- 0
  A[3, 1, year[t]] <- 0
  A[3, 2, year[t]] <- grow.na[t]
  A[3, 3, year[t]] <- phi.a[t]
}

# data
for(t in 1:time){
  y[1, t] ~ dpois(x[1, t])
  y[2, t] ~ dpois(x[2, t])
  y[3, t] ~ dpois(x[3, t])
} # t

# first latent process this can be continuous (not pois, maybe gamma (0 bound))
x[1, 1] ~ dpois(1) 
x[2, 1] ~ dpois(1) 
x[3, 1] ~ dpois(1) 

for(t in 1:(time-1)){
  mu[1:3, t] <- A[,,year.index[t,1]] %*% x[1:3,t]
  x[1:3, t+1] ~ dmnorm(mu[1:3, t], SIGMA)
} # t

}"
load.module("glm")
jags.model <- jags.model(textConnection(tick.met),
                         data = data,
                         n.adapt = 100000,
                         n.chains = 1)

variable.names <- c("phi.l.mu",
                    "phi.n.mu",
                    "phi.a.mu",
                    "grow.ln.mu",
                    'grow.na.mu',
                    "repro.mu",
                    "phi.l",
                    "phi.n",
                    "phi.a",
                    "grow.ln",
                    'grow.na',
                    "repro",
                    "SIGMA",
                    "beta.pl",
                    "beta.pn",
                    "beta.pa",
                    "beta.gln",
                    "beta.gna",
                    "beta.r",
                    "x")

xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file

jags.out <- coda.samples(jags.model,
                         n.iter = 250000,
                         variable.names = variable.names)

saveRDS(jags.out, 
        file = paste("/projectnb/dietzelab/fosterj/Previous_Runs_Tick/CaryTick_AnnualMet_", xx, ".rds",
                     sep = ""))





