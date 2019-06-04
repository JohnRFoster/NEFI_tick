library(rjags)

dat <- read.csv("tick_cleaned") # tick data
t <- subset(dat, Grid == "Green Control")
t <- t[,c("n_larvae", "n_nymphs", "n_adults")]
dat.t <- t(t) ## data for MCMC

gr <- read.csv("GreenTicks.csv", header = FALSE) # green control mice data
ch <- as.matrix(read.csv("GreenCaptHist.csv")) # green capture history
ch <- apply(ch, 2, as.numeric)

ks <- as.matrix(read.csv("KnownStatesGreen.csv")) # green capture history known states
ks <- apply(ks, 2, as.numeric)

min.alive <- apply(ks, 2, sum) # minimum number alive

### tick data
green <- dat[,c("Grid", "DATE", "n_larvae", "n_nymphs", "n_adults")]
green$Grid <- as.character(green$Grid)
green <- subset(green, Grid == "Green Control")
green$DATE <- as.Date(green$DATE)

### mouse data
mice <- read.csv("Cary_mouse.csv")
mice <- subset(mice, Grid == "Green Control")
mice$Full.Date.1 <- as.Date(mice$Full.Date.1)
mice$Full.Date.2 <- as.Date(mice$Full.Date.2)

# unique trap 
day.1 <- unique(mice$Full.Date.1) # 1st capture event in mouse series
day.1 <- as.character(day.1)
day.2 <- unique(mice$Full.Date.2) # 2nd capture event in mouse series
day.2 <- as.character(day.2)

day.mice <- c(rbind(day.1, day.2)) # unique sampling days: mice
day.mice <- as.Date(day.mice)
day.tick <- green$DATE # unique sampling days: tick

mna.date <- cbind(as.character(day.mice), min.alive) ## min alive each sampling date
row.names(mna.date) <- 1:nrow(mna.date) ## fix row names to numbers

## mice counts one year previous to tick sampling occasion
year.t1 <- vector()
for(t in 1:length(day.tick)){
  up <- day.tick[t] - 365 + 14        ## mice trapping date one year previous +/- 14 days
  down <- day.tick[t] - 365 - 14
  for(m in 1:length(day.mice)){
    if((day.mice[m] >= down) && (day.mice[m] <= up)){
      year.t1[t] <- as.character(day.mice[m])      
    } 
  }
}
## match sampling date with corrisponding MNA
mna.t1 <- data.frame(nrow = 72, ncol = 2)
for(t in 1:length(year.t1)){
  mna.t1[t,1] <- year.t1[t]
  if(is.na(year.t1[t])){
    mna.t1[t,2] <- NA
  }
  if(!is.na(year.t1[t])){
    xx <- which(year.t1[t] == mna.date[,1])
    mna.t1[t,2] <- mna.date[xx,2]
  }
}
colnames(mna.t1) <- c("Date_t1", "MNA_t1")

## mice counts two years previous to tick sampling occasion
year.t2 <- vector()
for(t in 1:length(day.tick)){
  up <- day.tick[t] - 365*2 + 14      ## mice trapping date two years previous +/- 14 days
  down <- day.tick[t] - 365*2 - 14
  for(m in 1:length(day.mice)){
    if((day.mice[m] >= down) && (day.mice[m] <= up)){
      year.t2[t] <- as.character(day.mice[m])      
    } 
  }
}

## match sampling date with corrisponding MNA
mna.t2 <- data.frame(nrow = 72, ncol = 2)
for(t in 1:length(year.t2)){
  mna.t2[t,1] <- year.t2[t]
  if(is.na(year.t2[t])){
    mna.t2[t,2] <- NA
  }
  if(!is.na(year.t2[t])){
    xx <- which(year.t2[t] == mna.date[,1])
    mna.t2[t,2] <- mna.date[xx,2]
  }
}
colnames(mna.t2) <- c("Date_t2", "MNA_t2")
lag <- cbind(as.character(day.tick), mna.t1, mna.t2) ## combinded data
mna.data <- lag[,c(3,5)]
mna.data <- apply(mna.data, 2, as.integer)
mna.t1.mis <- which(is.na(mna.data[,1]))
mna.t2.mis <- which(is.na(mna.data[,2]))

R <- diag(1, 3, 3)

data <- list(y = dat.t, 
             R = R,
             time = ncol(dat.t),
             mx = mna.data[,2],
             x.mis = mna.t2.mis,
             b0 = 0,   # regression beta means
             Vb = 0.001)   # regression beta precisions

inits <- function(){list(phi.l.mu = rnorm(1, 0, 2),
                         phi.n.mu = rnorm(1, 0, 2),
                         phi.a.mu = rnorm(1, 0, 2),
                         grow.ln.mu = rnorm(1, 0, 2),
                         grow.na.mu = rnorm(1, 0, 2),
                         repro.mu = rnorm(1, 0, 2),
                         beta.pl = rnorm(1, 0, 0.5),
                         beta.pn = rnorm(1, 0, 0.5),
                         beta.pa = rnorm(1, 0, 0.5),
                         beta.gln = rnorm(1, 0, 0.5),
                         beta.gna = rnorm(1, 0, 0.5),
                         beta.r = rnorm(1, 0, 0.5))}


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
for(i in x.mis){mx[i] ~ dunif(0, 145)}

for(t in 1:72) {
logit(phi.l[t]) <- phi.l.mu + alpha.pl[t] + beta.pl*mx[t]
logit(phi.n[t]) <- phi.n.mu + alpha.pn[t] + beta.pn*mx[t]
logit(phi.a[t]) <- phi.a.mu + alpha.pa[t] 
logit(grow.ln[t]) <- grow.ln.mu + alpha.gln[t] + beta.gln*mx[t]
logit(grow.na[t]) <- grow.na.mu + alpha.gna[t] + beta.gna*mx[t]
log(repro[t]) <- repro.mu + alpha.r[t]
alpha.pl[t] ~ dnorm(0, tau.pl)
alpha.pn[t] ~ dnorm(0, tau.pn)
alpha.pa[t] ~ dnorm(0, tau.pa)
alpha.gln[t] ~ dnorm(0, tau.gln)
alpha.gna[t] ~ dnorm(0, tau.gna)
alpha.r[t] ~ dnorm(0, tau.r)
}

# define transition matrix 
for(t in 1:72){
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
mu[1:3, t] <- A[,,t] %*% x[1:3,t]
x[1:3, t+1] ~ dmnorm(mu[1:3, t], SIGMA)
} # t

}"
load.module("glm")
jags.model <- jags.model(textConnection(tick.met),
                         data = data,
                         inits = inits,
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

saveRDS(jags.out,
        file = paste("/projectnb/dietzelab/fosterj/Previous_Runs_Tick/CaryTick_MiceTwo__", xx, ".rds",
                     sep = ""))
