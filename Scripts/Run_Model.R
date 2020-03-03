library(rjags)

#### model ####
source("Models/GDD_Null_K_All.R") 

n.adapt <- 50000
n.chains <- 3

## model options
met.variable <- NULL          # met driver on survival

## file path to output folder
out.folder <- "../FinalOut/Independent_Fits/GDDThreshold/K_all"
#out.folder <- "../FinalOut/HB_Partial_GDD/K_estimate"


### Individual Runs ###
site.run <- "Tea Control"      # site
out.name <- paste("NULL_K_est_multiProc", gsub(" ","",site.run),sep="_")
out.path <- file.path(out.folder,out.name)
file <- paste(out.path, ".RData", sep = "")

out <- run_model(site.run, n.adapt, n.chains, file)


### Hierarchical Runs ###
# out.name <- "GDDSwitch_Low_k"
# out <- run_model(n.adapt,n.chains,burnin,thin,n.iter)


# xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file
save(out, file = file)

print("-----------------------------------------------------------------")
print("                           END SCRIPT                            ")
print("-----------------------------------------------------------------")
print("                                                                 ")