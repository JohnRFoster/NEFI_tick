library(rjags)

t <- read.csv("tick_cleaned")
t <- subset(t, Grid == "Green Control")
t <- t[,c("n_larvae", "n_nymphs", "n_adults")]
dat.t <- t(t)

R <- diag(1, 3, 3)

data <- list(y = dat.t, 
             R = R,
             time = ncol(dat.t))

inits <- function(){list(x = dat.t,
                         phi.l = runif(1, 0, 1),
                         phi.n = runif(1, 0, 1),
                         phi.a = runif(1, 0, 1),
                         grow.ln = runif(1, 0, 1),
                         grow.na = runif(1, 0, 1))}

tick.null = " 
model {

# priors
phi.l ~ dunif(0, 1)       # larvae survival
phi.n ~ dunif(0, 1)       # nymph survival
phi.a ~ dunif(0, 1)       # adult survival
grow.ln ~ dunif(0, 1)     # larvae -> nymph transition (*truncated(phi.l))
grow.na ~ dunif(0, 1)     # nymph -> adult transition
repro ~ dgamma(1, 1)     # adult -> larvae transition (reproduction)
SIGMA ~ dwish(R, 3)       # variance matrix for mvn [3 x 3] default SIGMA will be *precision*
p.1 ~ dunif(0, 1)         # probability in larvae data model
p.2 ~ dunif(0, 1)         # probability in nymph data model
p.3 ~ dunif(0, 1)         # probability in adult data model

# define transition matrix - prior on whole matrix??
A[1, 1] <- phi.l
A[1, 2] <- 0
A[1, 3] <- repro
A[2, 1] <- grow.ln
A[2, 2] <- phi.n
A[2, 3] <- 0
A[3, 1] <- 0
A[3, 2] <- grow.na
A[3, 3] <- phi.a

# data
for(t in 1:time){
y[1, t] ~ dnegbin(p.1, x[1, t])
y[2, t] ~ dnegbin(p.2, x[2, t])
y[3, t] ~ dnegbin(p.3, x[3, t])
} # t

# first latent process this can be continuous (not pois, maybe gamma (0 bound))
x[1, 1] ~ dpois(1) 
x[2, 1] ~ dpois(1) 
x[3, 1] ~ dpois(1) 

for(t in 1:(time-1)){
mu[1:3, t] <- A %*% x[1:3,t]
x[1:3, t+1] ~ dmnorm(mu[1:3, t], SIGMA)
} # t

}"

jags.model <- jags.model(textConnection(tick.null),
                         data = data,
                         n.chains = 1)

variable.names <- c("phi.l",
                    "phi.n",
                    "phi.a",
                    "grow.ln", 
                    "grow.na",
                    "repro",
                    "p.1",
                    "p.2",
                    "p.3",
                    "SIGMA",
                    "x")

xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file

jags.out <- coda.samples(jags.model,
                         n.iter = 250000,
                         variable.names = variable.names)

saveRDS(jags.out, 
        file = paste("/projectnb/dietzelab/fosterj/Previous_Runs_Tick/CaryTickNull_NegBinomAll_", xx, ".rds",
                     sep = ""))





