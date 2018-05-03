library(rjags)


test1 <- readRDS("Cary_pop/CaryTest_Null_Out1.rds")
test2 <- readRDS("Cary_pop/CaryTest_Null_Out2.rds")
test3 <- readRDS("Cary_pop/CaryTest_Null_Out3.rds")


null.ls <- mcmc.list(c(test1,test2,test3))




