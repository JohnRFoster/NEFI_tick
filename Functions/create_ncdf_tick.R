library(ncdf4)

create_ncdf_tick <- function(ncfname, preds, params, start.date, data.assimilation,  
                             Nmc = 5000, 
                             ForecastProject.id = NULL, 
                             Forecast.id = NULL,
                             Forecast.issue.time = NULL){

  # ncfname = name of file to save
  # preds = matrix of state predictions (iterations x states)
  # params = matrix of parameters (iterations x parameters)
  # start.date = day forecast starts
  # data.assimilation = vector of 1s and 0s defining when we assimilate data
  # Nmc = number of samples to save (iterations)
  # ForecastProject.id = project identifier
  # Forecast.id = forecast identifier
  # Forecast.issue.time = date forecast was run
  
  # first let's thin samples to Nmc
  if(nrow(preds) > Nmc){
    thin <- round(seq(1, nrow(preds), length.out = Nmc))
    preds <- preds[thin,]
  }
  if(nrow(params) > Nmc){
    thin <- round(seq(1, nrow(params), length.out = Nmc))
    params <- params[thin,]
  }
  
  # need to split preds in to life stages
  larva.pred <- preds[,grep("x[1,", colnames(preds), fixed = TRUE)]
  nymph.pred <- preds[,grep("x[2,", colnames(preds), fixed = TRUE)]
  adult.pred <- preds[,grep("x[3,", colnames(preds), fixed = TRUE)]
  rm(preds)
  
  # number of days in forecast
  # remember the first day is the analysis day
  n.days <- ncol(larva.pred)
  
  # convert to date
  start.date <- ymd(start.date)
  
  # days forecasted
  time <- seq.Date(from = start.date, by = 1, length.out = n.days, format = "%Y")
  
  # Set dimensions
  ens <- 1:nrow(larva.pred)   # number of samples
  states <- 1                 # each life stage gets own dim
  timestep <- seq_along(time) # number days forecasted + initial condition
  n.params <- 1:ncol(params)  # number of parameters 
  max.param.char <- max(nchar(as.character(colnames(params)))) # max char length for param names

  sample.dim <- ncdim_def("ens", 
                      units = "mcmc_sample",
                      vals = ens, 
                      longname = 'posterior_sample') 
  state.dim <- ncdim_def("state", 
                        units = "individuals",
                        vals = states, 
                        longname = 'predicted_individuals') 
  time.dim <- ncdim_def("timestep", 
                       units = '1_day', 
                       longname = 'timestep',
                       vals = timestep)
  param.dim <- ncdim_def("parameters",
                         units = "",
                         vals = n.params,
                         longname = "parameter posteriors")
  da.dim <- ncdim_def("data_assimilation",
                      units = "boolean",
                      longname = "timestep when data is assimilated",
                      vals = timestep)
  dimnchar   <- ncdim_def("nchar",   "", 
                          1:nchar(as.character(time[1])), 
                          create_dimvar = FALSE)
  dimnchar.params   <- ncdim_def("nchar_params",   "",
                          1:max.param.char,
                          create_dimvar = FALSE)
  
  # Define variables
  fillvalue <- 1e32 # missing value
  
  def_list <- list()
  def_list[[1]] <- ncvar_def(name = "time",
                             units = "datetime",
                             dim = list(dimnchar, time.dim),
                             longname = "time",
                             prec = "char")
  def_list[[2]] <- ncvar_def(name =  "larvae",
                             units = "number of individuals",
                             dim = list(state.dim, sample.dim, time.dim),
                             missval = fillvalue,
                             longname = 'number_larvae_predicted',
                             prec = "single")
  def_list[[3]] <- ncvar_def(name =  "nymph",
                             units = "number of individuals",
                             dim = list(state.dim, sample.dim, time.dim),
                             missval = fillvalue,
                             longname = 'number_nymph_predicted',
                             prec = "single")
  def_list[[4]] <- ncvar_def(name =  "adult",
                             units = "number of individuals",
                             dim = list(state.dim, sample.dim, time.dim),
                             missval = fillvalue,
                             longname = 'number_adult_predicted',
                             prec = "single")
  def_list[[5]] <- ncvar_def(name =  "parameter_samples",
                             units = "parameter dependent",
                             dim = list(sample.dim, param.dim),
                             missval = fillvalue,
                             longname = 'parameters_from_mcmc',
                             prec = "single")
  def_list[[6]] <- ncvar_def(name = "parameter_names",
                             units = "name",
                             dim = list(dimnchar.params, param.dim),
                             longname = "parameter_names",
                             prec = "char")
  def_list[[7]] <- ncvar_def(name =  "data_assimilation",
                             units = "boolean",
                             dim = list(da.dim),
                             missval = fillvalue,
                             longname = '1 = data assimilation used in timestep',
                             prec = "single")

  # create nc file
  ncout <- nc_create(ncfname, def_list, force_v4 = TRUE)
  
  # put data into respective dimensions
  ncvar_put(ncout, def_list[[1]], time)                # dates in forecast
  ncvar_put(ncout, def_list[[2]], larva.pred)          # larva forecasts (first = analysis)
  ncvar_put(ncout, def_list[[3]], nymph.pred)          # nymph forecasts (first = analysis)  
  ncvar_put(ncout, def_list[[4]], adult.pred)          # adult forecasts (first = analysis)
  ncvar_put(ncout, def_list[[5]], params)              # parameter estimates
  ncvar_put(ncout, def_list[[6]], colnames(params))    # parameter names
  ncvar_put(ncout, def_list[[7]], data.assimilation)   # data assimilation
 
  # Global file metadata
  if (!is.null(ForecastProject.id)) {
    ncatt_put(ncout,
              0,
              "ForecastProject_id",
              as.character(ForecastProject.id),
              prec =  "text")
  }
  if (!is.null(Forecast.id)) {
    ncatt_put(ncout, 
              0, 
              "Forecast_id", 
              as.character(Forecast.id),
              prec =  "text")
  }
  if (!is.null(Forecast.issue.time)) {
    ncatt_put(ncout,
              0,
              "forecast_issue_time",
              as.character(Forecast.issue.time),
              prec =  "text")
  }
  nc_close(ncout)  
}






