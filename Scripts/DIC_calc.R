library(ecoforecastR)

#### model ####
model <- "Models/GDD_K_set.R" 
source(model) 


## Independent Models ##
sites <- c("Green Control", "Henry Control", "Tea Control")
xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file
site.run <- sites[xx]
cat("\n---", site.run, "---\n\n")
met.proc <- NULL            # met driver on survival


## HB Models ##
# met.proc <- c("max temp", "max rh", "vpd")            # met driver on survival
# cat("\n---", met.proc[xx], "---\n\n")

n.adapt <- 300000                 # adaptive iterations
n.chains <- 5                     # number of chains
n.iter <- 10000                  # number of iterations post burnin



# return <- run_model(met.proc[xx], n.adapt, n.chains) # for HB fits
return <- run_model(site.run, met.proc, n.adapt, n.chains) # for independent fits with met.proc
# return <- run_model(site.run, n.adapt, n.chains) # for independent fits 

cat("DIC for", model, "at", site.run, "\n")
dic.samples(model = return$j.model,
            n.iter = n.iter)