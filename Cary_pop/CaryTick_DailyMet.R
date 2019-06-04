library(rjags)

t <- read.csv("Cary_pop/tick_cleaned")
t <- subset(t, Grid == "Green Control")
df <- as.Date(t$DATE)

met <- read.csv("Cary_pop/Met_Cary") # daily weather data
met$DATE <- as.Date(met$DATE)

met <- met[which(met$DATE==df[1]):which(met$DATE==df[72]),]

df <- as.Date(t$DATE)
df <- diff(df)
df <- c(1, df)
dt.index <- cumsum(as.integer(df))

t <- t[,c("n_larvae", "n_nymphs", "n_adults")]
dat.t <- t(t)

precip <- met[,"TOT_PREC"]  # daily precip

temp <- met[,"MAX_TEMP"]   # daily max temperature
temp <- scale(temp, scale = FALSE)  # centered daily max temperature; range (-30.48,23.61) 
temp.mis <- which(is.na(temp))  # vector of indices for missing temp values

rh <- met[, "MAX_RH"] # daily max relative humidity
rh <- scale(rh, scale = FALSE)
rh.mis <- which(is.na(rh))  # vector of indices for missing rh values

R <- diag(1, 3, 3)

data <- list(y = dat.t, 
             R = R,
             time = ncol(dat.t),
             days = length(precip),
             temp = temp,
             temp.mis = temp.mis,
             precip = precip,
             rh = rh,
             rh.mis = rh.mis,
             dt.index = dt.index,
             b0 <- as.vector(c(0,0,0)),
             solve(diag(10000, 3)))

inits <- function(){list(phi.l.mu = rnorm(1, 0, 0.5),
                         phi.n.mu = rnorm(1, 0, 0.5),
                         phi.a.mu = rnorm(1, 0, 0.5),
                         grow.ln.mu = rnorm(1, 0, 0.5),
                         grow.na.mu = rnorm(1, 0, 0.5),
                         repro.mu = rnorm(1, 0, 0.5))}

tick.met = " 
model {

# priors
phi.l.mu ~ dnorm(0, 0.001)       # larvae survival
phi.n.mu ~ dnorm(0, 0.001)       # nymph survival
phi.a.mu ~ dnorm(0, 0.001)       # adult survival
grow.ln.mu ~ dnorm(0, 0.001)     # larvae -> nymph transition (*truncated(phi.l))
grow.na.mu ~ dnorm(0, 0.001)     # nymph -> adult transition
repro.mu ~ dnorm(0, 0.001) T(0,)       # adult -> larvae transition (reproduction)
SIGMA ~ dwish(R, 3)           # variance matrix for mvn [3 x 3] default SIGMA will be *precision*
beta.pl ~ dmnorm(b0,Vb)       # prior regression params larvae survival
beta.pn ~ dmnorm(b0,Vb)       # prior regression params nymph survival
beta.pa ~ dmnorm(b0,Vb)       # prior regression params adult survival
beta.gln ~ dmnorm(b0,Vb)      # prior regression params larvae -> nymph transition
beta.gna ~ dmnorm(b0,Vb)      # prior regression params nymph -> adult transition
beta.r ~ dmnorm(b0,Vb)        # prior regression params reproduction
tau.pl ~ dgamma(3, 2)
tau.pn ~ dgamma(3, 2)
tau.pa ~ dgamma(3, 2)
tau.gln ~ dgamma(3, 2)
tau.gna ~ dgamma(3, 2)
tau.r ~ dgamma(3, 2)
for(i in 1:temp.mis){temp[i] ~ dunif(-31, 24)}  # prior on missing temp
for(i in 1:rh.mis){rh[i] ~ dunif(-53, 6)}  # prior on missing rh

for(t in days) {
logit(phi.l[t]) <- phi.l.mu + beta.pl[1]*precip[t] + beta.pl[2]*temp[t] + beta.pl[3]*rh[t] + alpha.pl[t]
logit(phi.n[t]) <- phi.n.mu + beta.pn[1]*precip[t] + beta.pn[2]*temp[t] + beta.pn[3]*rh[t] + alpha.pn[t] 
logit(phi.a[t]) <- phi.a.mu + beta.pa[1]*precip[t] + beta.pa[2]*temp[t] + beta.pa[3]*rh[t] + alpha.pa[t] 
logit(grow.ln[t]) <- grow.ln.mu + beta.gln[1]*precip[t] + beta.gln[2]*temp[t] + beta.gln[3]*rh[t] + alpha.gln[t] 
logit(grow.na[t]) <- grow.na.mu + beta.gna[1]*precip[t] + beta.gna[2]*temp[t] + beta.gna[3]*rh[t] + alpha.gna[t]
log(repro[t]) <- repro.mu + beta.r[1]*precip[t] + beta.r[2]*temp[t] + beta.r[3]*rh[t] + alpha.r[t]
alpha.pl[t] ~ dnorm(0, tau.pl)
alpha.pn[t] ~ dnorm(0, tau.pn)
alpha.pa[t] ~ dnorm(0, tau.pa)
alpha.gln[t] ~ dnorm(0, tau.gln)
alpha.gna[t] ~ dnorm(0, tau.gna)
alpha.r[t] ~ dnorm(0, tau.r)
}

# define transition matrix 
for(t in days){
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
mu[1:3, t] <- A[,,yr_mon.index[t]] %*% x[1:3,t]
x[1:3, t+1] ~ dmnorm(mu[1:3, t], SIGMA)
} # t

}"

load.module("glm")

jags.model <- jags.model(textConnection(tick.met),
                         data = data,
                         inits = inits,
                         n.chains = 1)

variable.names <- c("phi.l",
                    "phi.n",
                    "phi.a",
                    "grow.ln", 
                    "grow.na",
                    "repro",
                    "p.1",
                    "SIGMA",
                    "temp.mis",
                    "rh.mis",
                    "x",
                    "beta.pl",
                    "beta.pn",
                    "beta.pa",
                    "beta.gln",
                    "beta.gna",
                    "beta.r")

#xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file

jags.out <- coda.samples(jags.model,
                         n.iter = 100,
                         variable.names = variable.names)

saveRDS(jags.out, 
        file = paste("/projectnb/dietzelab/fosterj/Previous_Runs_Tick/CaryTickNull_NegBinomLarv_", xx, ".rds",
                     sep = ""))





