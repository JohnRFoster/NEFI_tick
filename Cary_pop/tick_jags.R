



tick.model = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dnorm(x[i], tau_obs)
}

#### Process Model nymph
for(i in 2:n){
  Tnew[i, t] <- x[i, t-1] + beta.doy*doy[150, ] + beta.rh*rh[] + beta.precip*precip[] 

  x[i,t] ~ dnorm(Tnew[i, t], tau_add)
}



### Day of Year
for(t in 1:doy){
  doy[t] ~ dnorm(152, )
}

#### Priors
x[1] ~ dnorm(x_ic, tau_ic)
tau_obs ~ dgamma(a_obs, r_obs)
tau_add ~ dgamma(a_add, r_add)
}
"
### Green Control ###
data.green <- list(y = nymph.green, 
                   n = length(nymph.green), 
                   x_ic = 0.06, tau_ic = 1/0.06, 
                   a_obs = 1, r_obs = 0.01, 
                   a_add = 1, r_add = 0.01)