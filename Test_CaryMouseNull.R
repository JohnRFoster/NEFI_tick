
library(rjags)



j.model <- readRDS("Cary_Null_Model.rds")

jags.out <- coda.samples(model = j.model,
                         variable.names = c("lambda", "theta", "N"),
                         n.iter = 50000,
                         thin = 10)

saveRDS(jags.out, file = "Cary_Null_Out.rds")


