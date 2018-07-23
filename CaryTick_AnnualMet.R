library(rjags)
library(ecoforecastR)
t <- read.csv("tick_cleaned")

t <- subset(t, Grid == "Green Control")

ymd <- as.data.frame(strsplit(as.character(t$DATE),"-"))
year.index <- t(ymd[1,])
year.index <- apply(year.index, 2, as.integer)
year.index <- transform(year.index,id=as.numeric(factor(year.index[,1])))
year.index <- year.index[,2]

t <- t[,c("n_larvae", "n_nymphs", "n_adults")]
dat.t <- t(t)
dat.t <- apply(dat.t, 2, as.integer)

### average annual
year <- 1995:2005
met <- read.csv("Met_Cary")
met <- met[, c("DATE","MAX_TEMP","MAX_RH","TOT_PREC")]
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

temp <- scale(temp, scale = FALSE) ## centered mean annual temp
temp <- temp[,1]

rh <- scale(rh, scale = FALSE) ## centered mean annual rh
rh <- rh[,1]

precip <- scale(precip, scale = FALSE)
precip <- precip[,1]

R <- diag(1, 3, 3)


year <- 1:11 
data <- list(y = dat.t, 
             R = R,
             time = ncol(dat.t),
             year.index = year.index,
             temp = temp,
             precip = precip,
             rh = rh,
             b0 = as.vector(c(0,0,0)),   # regression beta means
             Vb = solve(diag(10000,3)))   # regression beta precisions

noise <- function(dat, sd){
  for(c in 1:ncol(dat)){
    for(r in 1:nrow(dat)){
      dat[r, c] <- round(rnorm(n = 1, mean = dat[r, c], sd = sd))
      if(dat[r, c] < 0){dat[r, c] <- 0}
    }
  }
  return(dat)
}

x <- noise(dat.t, 0.5)

inits <- function(){list(phi.l.mu = rnorm(1, 3.564, 0.5),
                         phi.n.mu = rnorm(1, 1.61, 0.2),
                         phi.a.mu = rnorm(1, -0.13, 0.1),
                         grow.ln.mu = rnorm(1, -1.27, 0.1),
                         grow.na.mu = rnorm(1, -0.069, 0.15),
                         repro.mu = rnorm(1, -0.95, 0.15))}
# beta.pl = rnorm(1, 0, 0.01),
# beta.pn = rnorm(1, 0, 0.01),
# beta.pa = rnorm(1, 0, 0.01),
# beta.gln = rnorm(1, 0, 0.01),
# beta.gna = rnorm(1, 0, 0.01),
# beta.r = rnorm(1, 0, 0.01))}
#x = x)}

tick.met = " 
model {

# priors
phi.l.mu ~ dnorm(0, 0.3)       # larvae survival
phi.n.mu ~ dnorm(0, 0.3)       # nymph survival
phi.a.mu ~ dnorm(0, 0.3)       # adult survival
grow.ln.mu ~ dnorm(0, 0.3)     # larvae -> nymph transition (*truncated(phi.l))
grow.na.mu ~ dnorm(0, 0.3)     # nymph -> adult transition
repro.mu ~ dnorm(0, 0.3)       # adult -> larvae transition (reproduction)
SIGMA ~ dwish(R, 4)           # variance matrix for mvn [3 x 3] default SIGMA will be *precision*
beta.pl ~ dmnorm(b0, Vb)       # prior regression params larvae survival
beta.pn ~ dmnorm(b0, Vb)       # prior regression params nymph survival
beta.pa ~ dmnorm(b0, Vb)       # prior regression params adult survival
beta.gln ~ dmnorm(b0, Vb)      # prior regression params larvae -> nymph transition
beta.gna ~ dmnorm(b0, Vb)      # prior regression params nymph -> adult transition
beta.r ~ dmnorm(b0, Vb)        # prior regression params reproduction
tau.pl ~ dgamma(3, 2)
tau.pn ~ dgamma(3, 2)
tau.pa ~ dgamma(3, 2)
tau.gln ~ dgamma(3, 2)
tau.gna ~ dgamma(3, 2)
tau.r ~ dgamma(3, 2)


for(t in 1:11) {
logit(phi.l[t]) <- phi.l.mu + alpha.pl[t] + beta.pl[1]*rh[t] + beta.pl[2]*precip[t] + beta.pl[3]*temp[t]
logit(phi.n[t]) <- phi.n.mu + alpha.pn[t] + beta.pn[1]*rh[t] + beta.pn[2]*precip[t] + beta.pn[3]*temp[t]
logit(phi.a[t]) <- phi.a.mu + alpha.pa[t] + beta.pa[1]*rh[t] + beta.pa[2]*precip[t] + beta.pa[3]*temp[t]
logit(grow.ln[t]) <- grow.ln.mu + alpha.gln[t] + beta.gln[1]*rh[t] + beta.gna[2]*precip[t] + beta.gln[3]*temp[t] 
logit(grow.na[t]) <- grow.na.mu + alpha.gna[t] + beta.gna[1]*rh[t] + beta.gna[2]*precip[t] + beta.gna[3]*temp[t]
log(repro[t]) <- repro.mu + alpha.r[t] + beta.r[1]*rh[t] + beta.r[2]*precip[t] + beta.r[3]*temp[t]
alpha.pl[t] ~ dnorm(0, tau.pl)
alpha.pn[t] ~ dnorm(0, tau.pn)
alpha.pa[t] ~ dnorm(0, tau.pa)
alpha.gln[t] ~ dnorm(0, tau.gln)
alpha.gna[t] ~ dnorm(0, tau.gna)
alpha.r[t] ~ dnorm(0, tau.r)
}

# define transition matrix 
for(t in 1:11){
A[1, 1, t] <- phi.l[t]
A[1, 2, t] <- 0
A[1, 3, t] <- repro[t]
A[2, 1, t] <- grow.ln[t]
A[2, 2, t] <- phi.n[t]
A[2, 3, t] <- 0
A[3, 1, t] <- 0
A[3, 2, t] <- grow.na[t]
A[3, 3, t] <- phi.a[t]
}

# data
for(t in 1:time){
y[1, t] ~ dpois(x[1, t])
y[2, t] ~ dpois(x[2, t])
y[3, t] ~ dpois(x[3, t])
} # t



# first latent process this can be continuous (not pois, maybe gamma (0 bound))
x[1, 1] ~ dpois(2)
x[2, 1] ~ dpois(1)
x[3, 1] ~ dpois(9)

for(t in 1:(time-1)){
mu[1:3, t] <- A[,,year.index[t]] %*% x[1:3,t]
x[1:3, t+1] ~ dmnorm(mu[1:3, t], SIGMA)
} # t

}"
load.module("glm")
jags.model <- jags.model(textConnection(tick.met),
                         data = data,
                         inits = inits,
                         n.adapt = 100000,
                         n.chains = 1)

#dic.samples(jags.model, n.iter = 250000, n.adapt = 100000)

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
                    "x",
                    "alpha.pl",
                    "alpha.pn",
                    "alpha.pa",
                    "alpha.gln",
                    "alpha.gna",
                    "alpha.r",
                    "tau.pl",
                    "tau.pa",
                    "tau.pn",
                    "tau.gln",
                    "tau.gna",
                    "tau.r")



xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file

jags.out <- coda.samples(jags.model,
                         n.iter = 250000,
                         variable.names = variable.names)

#jags.out <- window(jags.out, start = 400)
## split output

mfit <- as.matrix(jags.out, chains = TRUE)


phi.l.mu <- mfit[, grep("phi.l.mu",colnames(mfit))]
print("phi.l.mu")
quantile(phi.l.mu, c(0.25, 0.5, 0.75))

phi.n.mu <- mfit[, grep("phi.a.mu",colnames(mfit))]
print("phi.n.mu")
quantile(phi.n.mu, c(0.25, 0.5, 0.75))

phi.a.mu <- mfit[, grep("phi.n.mu",colnames(mfit))]
print("phi.a.mu")
quantile(phi.a.mu, c(0.25, 0.5, 0.75))

grow.ln.mu <- mfit[, grep("grow.ln.mu",colnames(mfit))]
print("grow.ln.mu")
quantile(grow.ln.mu, c(0.25, 0.5, 0.75))

grow.na.mu <- mfit[, grep("grow.na.mu",colnames(mfit))]
print("grow.na.mu")
quantile(grow.na.mu, c(0.25, 0.5, 0.75))

repro.mu <- mfit[, grep("repro.mu",colnames(mfit))]
print("repro.mu")
quantile(repro.mu, c(0.25, 0.5, 0.75))

saveRDS(jags.out,
        file = paste("/projectnb/dietzelab/fosterj/Previous_Runs_Tick/CaryTick_AnnualMet_", xx, ".rds",
                     sep = ""))





