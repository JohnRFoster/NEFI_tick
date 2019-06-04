print("                                                                 ")
print("-----------------------------------------------------------------")
print("                         Start SCRIPT                            ")
print("-----------------------------------------------------------------")


library(ecoforecastR)

source("Functions/cary_tick_met_JAGS.R")
source("Functions/site_data_met.R")
source("FinalOut/Models/Temp_MoltingGDD_Independent.R")

## model options
site.run <- "Green Control"
met.variable <- "temp"
n.adapt <- 10000                  # adpative iterations
n.chains <- 1                     # number of chains
burnin <- 25000                   # burnin iterations
n.iter <- 10000                  # number of iterations post burnin
iter2save <- 5000                 # number of iterations to save 
thin <- round(n.iter/iter2save)   # thinning interval

## file path to output folder
out.folder <- "FinalOut/Independent_Fits/WithMolting"
out.name <- paste("Temp_MoltingeEXP", gsub(" ","",site.run),sep="_")
out.path <- file.path(out.folder,out.name)

# compile and run model
out <- run_model(site.run,met.variable,n.adapt,n.chains,burnin,thin,n.iter)

xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file
save(out, file = paste(out.path, xx, ".RData", sep = ""))

print("-----------------------------------------------------------------")
print("                           END SCRIPT                            ")
print("-----------------------------------------------------------------")
print("                                                                 ")