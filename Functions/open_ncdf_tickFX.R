# this script has functions for openning ncdf files generated from the tick forecasts


# for data assimilation in the forecast
get_fx_da <- function(dir, day.2.grab){
  
  files.forecast <- list.files(dir)
  
  # assimilate forecast from climate ensembles
  if(any(grepl("clim", files.forecast))){
    files.forecast <- files.forecast[grep("clim", files.forecast)]
    nc <- nc_open(file.path(dir, files.forecast[1])) # read nc file
    nrow <- nrow(ncvar_get(nc, "larvae")) # for dimensions
    params.nc <- matrix(0, 1, length(params.nc.names))
    colnames(params.nc) <- ncvar_get(nc, "parameter_names") # don't change
    nc.larva <- nc.nymph <- nc.adult <- matrix(NA, nrow, length(files.forecast))
    for(cf in seq_along(files.forecast)){
      nc <- nc_open(file.path(dir, files.forecast[cf])) # read nc file
      f.days <- ncvar_get(nc, "time")     # dates for each forecast
      f.2.grab <- which(day.2.grab == f.days) # the day we need
      nc.larva[,cf] <- ncvar_get(nc, "larvae")[,f.2.grab] # larva forecasts
      nc.nymph[,cf] <- ncvar_get(nc, "nymph")[,f.2.grab]  # nymph forecasts
      nc.adult[,cf] <- ncvar_get(nc, "adult")[,f.2.grab]  # adult forecasts
      params.ens <- ncvar_get(nc, "parameter_samples")
      colnames(params.ens) <- ncvar_get(nc, "parameter_names") # don't change
      params.nc <- rbind(params.nc, params.ens)
    }
    params.nc <- params.nc[-1,] # remove first row from initialization above
  } else {
    nc <- nc_open(file.path(dir, files.forecast[1])) # read nc file
    f.days <- ncvar_get(nc, "time")     # dates for each forecast
    f.2.grab <- which(day.2.grab == f.days) # the day we need
    nc.larva <- ncvar_get(nc, "larvae")[,f.2.grab] # larva forecasts
    nc.nymph <- ncvar_get(nc, "nymph")[,f.2.grab]  # nymph forecasts
    nc.adult <- ncvar_get(nc, "adult")[,f.2.grab]  # adult forecasts
    params.nc <- ncvar_get(nc, "parameter_samples")
    colnames(params.nc) <- ncvar_get(nc, "parameter_names") # don't change
  }
  
  # pull out median predictions for DA
  larva.m <- median(nc.larva)
  nymph.m <- median(nc.nymph)
  adult.m <- median(nc.adult)
  
  
  return(list(
    ic = c(larva.m,
           nymph.m,
           adult.m),
    params.nc = params.nc
  ))

}
