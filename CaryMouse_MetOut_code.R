library(rjags)
library(runjags)
library(ecoforecastR)

chain1 <- list()
for(i in 1:10){
  chain1[i] <- readRDS(file = paste("CaryMax_Recruit_Out_3.", i, ".rds", sep = ""))
}
c1 <- combine.mcmc(chain)1
print("chain 1")
chain2 <- list()
for(i in 1:10){
  chain2[i] <- readRDS(file = paste("CaryMax_Recruit_Out_2.", i, ".rds", sep = ""))
}
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

effect.size.Recruit <- effectiveSize(out$params)
print("effective sample size")

plot.Recruit <- plot(out$params)
print("plot")

trace.Recruit <- traceplot(out$params)
print("trace plot")

## grab params of interest
lambda.mean.Recruit <- mfit[, grep("lambda.mean", colnames(mfit))] 
print("grep lambda.mean")
N.Recruit <- mfit[, grep("N",colnames(mfit))]
print("grep N")
lambda.Recruit <- mfit[, grep("lambda",colnames(mfit))]
print("grep lambda")
n.mean.Recruit <- apply(N, 2, mean) # calculated abundance
print("calculate n.mean")

x <- read.csv("KnownStatesGreen.csv") # known states for each individual
x <- apply(x, 2, as.numeric)
n.caught <- apply(x, 2, sum) # minimum number alive
print("n.caught")

nsamp <- 5000
samp <- sample.int(nrow(mfit),nsamp)

ci.N.Recruit <- apply(N[samp,],2,quantile,c(0.025,0.5,0.975))
print("ci.N")

save(lambda.mean.Recruit,
     N.Recruit,
     lambda.Recruit,
     n.mean.Recruit,
     #n.caught,
     effect.size.Recruit,
     plot.Recruit,
     trace.Recruit,
     ci.N.Recruit,
     file = "/projectnb/dietzelab/fosterj/CaryMouseRecruit.RData")

prin("END FILE")