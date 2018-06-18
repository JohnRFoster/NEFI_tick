library(rjags)

dat.t <- read.csv("GreenTicks.csv", header = FALSE)

R <- diag(1, 3, 3)

data <- list(y = dat.t, 
             R = R,
             time = ncol(dat.t))

data$b0 <- as.vector(c(0, 0, 0))           # regression beta means
data$Vb <- solve(diag(10000, 3))           # regression beta precisions


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
repro ~ dgamma(1, 1)      # adult -> larvae transition (reproduction)
SIGMA ~ dwish(R, 3)       # variance matrix for mvn [3 x 3] default SIGMA will be *precision*

beta.phi.l ~ dmnorm(b0,Vb)                    # prior regression params
beta.phi.n ~ dmnorm(b0,Vb)                    # prior regression params
beta.phi.a ~ dmnorm(b0,Vb)                    # prior regression params
beta.grow.ln ~ dmnorm(b0,Vb)                  # prior regression params
beta.grow.na ~ dmnorm(b0,Vb)                  # prior regression params
beta.repro ~ dmnorm(b0,Vb)                    # prior regression params

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

log(phi.l) <- beta.phi.l[1] + beta.phi.l[2]*mice.t1[t-1] + beta.phi.l[3]*mice.t2 
log(phi.n) <- beta.phi.n[1] + beta.phi.n[2]*mice.t1[t-1] + beta.phi.n[3]*mice.t2 
log(phi.a) <- beta.phi.a[1] + beta.phi.a[2]*mice.t1[t-1] + beta.phi.a[3]*mice.t2 
log(grow.ln) <- beta.grow.ln[1] + beta.grow.ln[2]*mice.t1[t-1] + beta.grow.ln[3]*mice.t2 
log(grow.na) <- beta.grow.na[1] + beta.grow.na[2]*mice.t1[t-1] + beta.grow.na[3]*mice.t2 
log(repro) <- beta.repro[1] + beta.repro[2]*mice.t1[t-1] + beta.repro[3]*mice.t2 

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
mu[1:3, t] <- A %*% x[1:3,t]
x[1:3, t+1] ~ dmnorm(mu[1:3, t], SIGMA)
} # t

}"

jags.model <- jags.model(textConnection(tick.null),
                         data = data,
                         n.chains = 3)

variable.names <- c("x", 
                    "phi.l",
                    "phi.n",
                    "phi.a",
                    "grow.ln", 
                    "grow.na",
                    "repro",
                    "SIGMA")

jags.out <- coda.samples(jags.model,
                         n.iter = 1000000,
                         variable.names = variable.names)

saveRDS(jags.out, file = "CaryTickNull.rds")

