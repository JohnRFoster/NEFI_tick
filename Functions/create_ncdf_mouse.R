library(ncdf4)

create_ncdf_mouse <- function(ncfname, n, start.date, n.days, Nmc, data_assimilation,
                              ForecastProject_id, Forecast_id, forecast_issue_time){
  
  # convert to date
  start.date <- ymd(start.date)
  
  # days forecasted
  time <- seq.Date(from = start.date, by = 1, length.out = n.days, format = "%Y")
  
  # did we assimilate data
  data_assimilation <- rep(data_assimilation, length(time))
  
  # Set dimensions
  ens <- as.integer(seq(1, Nmc, 1)) # number of ensembles
  states <- as.integer(1:dim(n)[2]) # individuals
  timestep <- as.integer(1:n.days)
  
  ensdim <- ncdim_def("ens", 
                      units = "",
                      vals = ens, 
                      longname = 'ensemble_member') 
  statedim <- ncdim_def("state", 
                        units = "individuals",
                        vals = states, 
                        longname = 'predicted_individuals') 
  timedim <- ncdim_def("timestep", 
                       units = '1 day', 
                       longname = 'timestep',
                       vals = timestep)
  
  dimnchar   <- ncdim_def("nchar",   "", 
                          1:nchar(as.character(time[1])), 
                          create_dimvar = FALSE)
  
  #Define variables
  fillvalue <- 1e32 # missing value
  
  def_list <- list()
  def_list[[1]] <- ncvar_def(name = "time",
                             units = "datetime",
                             dim = list(dimnchar, timedim),
                             longname = "time",
                             prec ="char")
  def_list[[2]] <- ncvar_def(name =  "latent_state",
                             units = "number of individuals",
                             dim = list(statedim, ensdim, timedim),
                             missval = fillvalue,
                             longname = 'number_mice_alive',
                             prec ="single")
  def_list[[3]] <- ncvar_def(name =  "predicted_observed",
                             units = "number of individuals",
                             dim = list(statedim, ensdim, timedim),
                             missval = fillvalue,
                             longname = 'number_mice_observed_given_number_alive',
                             prec ="single")
  def_list[[4]] <- ncvar_def(name =  "data_assimilation",
                             units = "logical",
                             dim = list(timedim),
                             missval = fillvalue,
                             longname = '1 = data assimilation used in timestep',
                             prec = "single")
  
  
  ncout <- nc_create(ncfname, def_list, force_v4 = TRUE)
  
  
  
  ncvar_put(ncout,def_list[[1]] , time)
  ncvar_put(ncout,def_list[[2]] , n[1, , , ])
  ncvar_put(ncout,def_list[[3]] , n[2, , , ])
  ncvar_put(ncout,def_list[[4]] , data_assimilation)
  
  #Global file metadata
  ncatt_put(ncout,0,"ForecastProject_id", as.character(ForecastProject_id), 
            prec =  "text")
  ncatt_put(ncout,0,"Forecast_id",as.character(Forecast_id), 
            prec =  "text")
  ncatt_put(ncout,0,"forecast_issue_time",as.character(forecast_issue_time), 
            prec =  "text")
  nc_close(ncout)  
}
  





