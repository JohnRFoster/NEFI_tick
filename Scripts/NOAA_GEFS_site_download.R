# # need to set this for CRON, will default to geo home instead of projectnb
setwd("/projectnb/dietzelab/fosterj/NEFI_tick")

# devtools::install_github("rqthomas/noaaGEFSpoint")
library(rNOMADS)
# Sys.setenv(PATH = "/projectnb/dietzelab/fosterj/grib2/wgrib2/")
library(rgdal)
library(noaaGEFSpoint)
library(tidyverse)
library(lubridate)
library(parallel)

print(sessionInfo())
cat("\n\n")
cat("GEFS Specifications:\n")


output.directory <- "/projectnb/dietzelab/fosterj/Data"
cat("Output directory", output.directory, "\n")

# site_file <- system.file("extdata", "noaa_download_site_list.csv", package = "noaaGEFSpoint")
# neon_sites <- read.csv(site_file) # read site file .csv
# colnames(neon_sites)[1] <- "site_name" # rename
# cary.info <- data.frame(site_name = "Cary Institute", # add Cary to neon_sites df
#                         site_id = "CARY",
#                         latitude = 41.7851,
#                         longitude = -73.7338)
# neon_sites <- bind_rows(neon_sites, cary.info)
# 
# # subset to tick (Ixodes and Amblyomma) prominent sites 
# neon_sites <- neon_sites %>% 
#   filter(site_id %in% c("HARV", "CARY", "BLAN", "ORNL", "SCBI", 
#                         "SERC", "KONZ", "OSBS", "TALL", "UKFS"))
# cat(nrow(neon_sites), "sites:", neon_sites$site_id, "\n")

# set up parallel, default is false
n.cores <- 6 # for testing
# n.cores <- as.numeric(Sys.getenv("NSLOTS")) # read number of parallel slots (cores) set in bash script
cat("Number of cores for download:", n.cores, "\n")

run.par <- FALSE
if(n.cores > 1) run.par <- TRUE
cat("Running in parallel:", run.par, "\n")

# method for download, default is point
method <- "point"
# if(nrow(neon_sites) > 3 | run.par) method <- "grid"
cat("Download method:", method, "\n")

overwrite <- FALSE
downscale <- FALSE
forecast.time <- "all"
forecast.date <- "all"

cat("Downscale:", downscale, "\n")
cat("Overwrite:", overwrite, "\n")
cat("Forecast time:", forecast.time, "\n")
cat("Forecast date:", forecast.date, "\n")


# cat("Running noaa_gefs_download_downscale() with method = point\n")
# noaaGEFSpoint::noaa_gefs_download_downscale(site_list = "CARY",
#                                             lat_list = 41.7851,
#                                             lon_list = -73.7338,
#                                             output_directory = file.path(output.directory, "GEFSpoint"),
#                                             forecast_time = forecast.time,
#                                             forecast_date = forecast.date,
#                                             downscale = downscale,
#                                             run_parallel = run.par,
#                                             num_cores = n.cores,
#                                             method = "point",
#                                             overwrite = overwrite)

cat("Running noaa_gefs_download_downscale() with method = grid, run_parallel = FALSE\n")
noaaGEFSpoint::noaa_gefs_download_downscale(site_list = "CARY",
                                            lat_list = 41.7851,
                                            lon_list = -73.7338,
                                            output_directory = file.path(output.directory, "GEFSgrid"), 
                                            forecast_time = forecast.time,
                                            forecast_date = forecast.date,
                                            downscale = downscale, 
                                            run_parallel = FALSE, 
                                            num_cores = 1, 
                                            method = "grid",
                                            overwrite = overwrite)

print(warnings())
