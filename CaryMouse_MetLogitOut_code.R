library(rjags)
library(runjags)
library(ecoforecastR)

chain1 <- list()
for(i in 1:8){
  chain1[i] <- readRDS(file = paste("Cary_DailyMetLogit_Out_1.", i, ".rds", sep = ""))
}
c1 <- combine.mcmc(chain1)
print("chain 1")
chain2 <- list()
for(i in 1:8){
  chain2[i] <- readRDS(file = paste("Cary_DailyMetLogit_Out_2.", i, ".rds", sep = ""))
}
c2 <- combine.mcmc(chain2)
print("chain 2")

jags.out <- mcmc.list(c1,c2)
print("jags.out")

## split output
out <- list(params = NULL, predict = NULL)
mfit <- as.matrix(jags.out, chains = TRUE)
pred.cols <- grep("N[", colnames(mfit), fixed = TRUE)
chain.col <- which(colnames(mfit) == "CHAIN")

out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
print("out$predict")
out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
print("out$params")

# GBR <- gelman.diag(out$params)
# print("GBR")
# 
# burnin <- GBR$last.iter[tail(which(any(GBR$shrink[,,2] > 1.1)),1)+1]
# cat("Burnin: ", burnin)
# 
# if(length(burnin) == 0) burnin = 1
# 
# ## remove burn-in
# jags.burn <- window(jags.out,start = burnin)
# print("Remove burnin")
# 
# jags.summary <- summary(jags.burn)
# print("summary output")
# 
# mfit <- as.matrix(jags.burn)
# print("jags matrix post burnin")

effect.size <- effectiveSize(out$params)
print("effective sample size")

plot <- plot(out$params)
print("plot")

trace <- traceplot(out$params)
print("trace plot")

## grab params of interest
lambda.mean <- mfit[, grep("lambda.mean", colnames(mfit))] 
print("grep lambda.mean")
N <- mfit[, grep("N",colnames(mfit))]
print("grep N")
lambda <- mfit[, grep("lambda",colnames(mfit))]
print("grep lambda")
n.mean <- apply(N, 2, mean) # calculated abundance
print("calculate n.mean")

x <- read.csv("KnownStatesGreen.csv") # known states for each individual
x <- apply(x, 2, as.numeric)
n.caught <- apply(x, 2, sum) # minimum number alive
print("n.caught")

ci.N <- apply(N[samp,],2,quantile,c(0.025,0.5,0.975))
print("ci.N")

save(lambda.mean,
     N,
     lambda,
     n.mean,
     n.caught,
     effect.size,
     plot,
     trace,
     ci.N,
     file = "/projectnb/dietzelab/fosterj/CaryMouseMetLogit.RData")

# nsamp <- 5000
# samp <- sample.int(nrow(mfit),nsamp)
# dim <- c(nsamp, nrow(ch), ncol(ch))
# x.precip <- 1:10                ## sequence of x values we're going to
# npred <- length(xpred)              ##      make predictions for
# ypred <- matrix(0.0,nrow=nsamp,ncol=npred)  ## storage for predictive interval
# ycred <- matrix(0.0,nrow=nsamp,ncol=npred)  ## storage for credible interval


