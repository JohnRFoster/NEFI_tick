print("                                                                 ")
print("-----------------------------------------------------------------")
print("                         Start SCRIPT                            ")
print("-----------------------------------------------------------------")


library(ecoforecastR)

#### model ####
source("Models/Temp_GDDThreshold_Independent.R") 

## model options
site.run <- "Green Control"      # site
met.variable <- NULL          # met driver on survival
n.adapt <- 150000                  # adaptive iterations
n.chains <- 1                     # number of chains
burnin <- 50000                   # burnin iterations
n.iter <- 750000                  # number of iterations post burnin
iter2save <- 10000                 # number of iterations to save 
thin <- round(n.iter/iter2save)   # thinning interval

## file path to output folder
out.folder <- "../FinalOut/Independent_Fits/GDDThreshold/LowThreshOnly"
out.name <- paste("Temp_GDDSwitch_gammaBeta33", gsub(" ","",site.run),sep="_")
out.path <- file.path(out.folder,out.name)

# compile and run model
out <- run_model(site.run,n.adapt,n.chains,burnin,thin,n.iter)

xx <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number to paste into output file
save(out, file = paste(out.path, xx, ".RData", sep = ""))

print("-----------------------------------------------------------------")
print("                           END SCRIPT                            ")
print("-----------------------------------------------------------------")
print("                                                                 ")