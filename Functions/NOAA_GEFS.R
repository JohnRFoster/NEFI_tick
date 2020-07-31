

NOAA_GEFS <- function(outfolder, lat.in, lon.in, sitename, 
                      start_date = Sys.time(), 
                      end_date = (as.POSIXct(start_date, tz="UTC") + lubridate::days(16)),
                      time.zone = tz,
                      overwrite = FALSE, verbose = FALSE, ...) {
  
  start_date <- as.POSIXct(start_date, tz = "UTC")
  end_date <- as.POSIXct(end_date, tz = "UTC")
  
  #It takes about 2 hours for NOAA GEFS weather data to be posted.  Therefore, if a request is made within that 2 hour window,
  #we instead want to adjust the start time to the previous forecast, which is the most recent one avaliable.  (For example, if
  #there's a request at 7:00 a.m., the data isn't up yet, so the function grabs the data at midnight instead.)
  if (abs(as.numeric(Sys.time() - start_date, units="hours")) <= 2) {
    start_date = start_date - lubridate::hours(2)
    end_date = end_date - lubridate::hours(2)
  }
  
  #Date/time error checking - Checks to see if the start date is before the end date
  if (start_date > end_date) {
    PEcAn.logger::logger.severe("Invalid dates: end date occurs before start date")
  } else if (as.numeric(end_date - start_date, units="hours") < 6) { #Done separately to produce a more helpful error message.
    PEcAn.logger::logger.severe("Times not far enough appart for a forecast to fall between them.  Forecasts occur every six hours; make sure start 
                                and end dates are at least 6 hours appart.")
  }
  
  #Set the end forecast date (default is the full 16 days)
  if (end_date > start_date + lubridate::days(16)) {
    end_date = start_date + lubridate::days(16)
    PEcAn.logger::logger.info(paste0("Updated end date is ", end_date))
  }
  
  #Round the starting date/time down to the previous block of 6 hours.  Adjust the time frame to match.
  forecast_hour = (lubridate::hour(start_date) %/% 6) * 6 #Integer division by 6 followed by re-multiplication acts like a "floor function" for multiples of 6
  increments = as.integer(as.numeric(end_date - start_date, units = "hours") / 6) #Calculating the number of forecasts between start and end dates.
  increments = increments + ((lubridate::hour(end_date) - lubridate::hour(start_date)) %/% 6) #These calculations are required to use the rnoaa package.
  
  end_hour = sprintf("%04d", ((forecast_hour + (increments * 6)) %% 24) * 100)  #Calculating the starting hour as a string, which is required type to access the 
  #data via the rnoaa package
  forecast_hour = sprintf("%04d", forecast_hour * 100)  #Having the end date as a string is useful later, too.
  
  #Recreate the adjusted start and end dates.
  start_date = as.POSIXct(paste0(lubridate::year(start_date), "-", lubridate::month(start_date), "-", lubridate::day(start_date), " ", 
                                 substring(forecast_hour, 1,2), ":00:00"), tz="UTC")
  
  end_date = start_date + lubridate::hours(increments * 6)
  
  
  #Bounds date checking
  #NOAA's GEFS database maintains a rolling 12 days of forecast data for access through this function.
  #We do want Sys.Date() here - NOAA makes data unavaliable days at a time, not forecasts at a time.
  NOAA_GEFS_Start_Date = as.POSIXct(Sys.Date(), tz="UTC") - lubridate::days(11)  #Subtracting 11 days is correct, not 12.
  
  #Check to see if start_date is valid. This must be done after date adjustment.
  if (as.POSIXct(Sys.time(), tz="UTC") < start_date || start_date < NOAA_GEFS_Start_Date) {
    PEcAn.logger::logger.severe(sprintf('Start date (%s) exceeds the NOAA GEFS range (%s to %s).',
                                        start_date,
                                        NOAA_GEFS_Start_Date, Sys.Date()))
  }
  
  if (lubridate::hour(start_date) > 23) {
    PEcAn.logger::logger.severe(sprintf("Start time %s is not a valid time", lubridate::hour(start_date)))
  }
  
  if (lubridate::hour(end_date) > 23) {  #Done separately from the previous if statement in order to have more specific error messages.
    PEcAn.logger::logger.severe(sprintf("End time %s is not a valid time", lubridate::hour(end_date)))
  }
  #End date/time error checking
  
  # End date/time error checking
  
  #################################################
  # NOAA variable downloading
  # Uses the rnoaa package to download data
  
  # We want data for each of the following variables. 
  noaa_var_names = c("Temperature_height_above_ground_ens",
                     # "foobar",
                     "Relative_humidity_height_above_ground_ens",
                     "Total_precipitation_surface_6_Hour_Accumulation_ens")
  
  cf_var_names = noaa_var_names
  
  #Downloading the data here.  It is stored in a matrix, where columns represent time in intervals of 6 hours, and rows represent
  #each ensemble member.  Each variable getxs its own matrix, which is stored in the list noaa_data.
  
  noaa_data <- list()
  for (i in 1:length(noaa_var_names)) {
    print(noaa_var_names[i])
    noaa.mat <- rnoaa::gefs(noaa_var_names[i], 
                                 lat.in, lon.in, 
                                 raw = TRUE, 
                                 time_idx = seq_len(increments), 
                                 forecast_time = forecast_hour, 
                                 date = format(start_date, "%Y%m%d"))$data
    noaa_data[[i]] <- noaa.mat
  }
  
  #Fills in data with NaNs if there happens to be missing columns.
  for (i in 1:length(noaa_var_names)) {
    if (!is.null(ncol(noaa_data[[i]]))) { # Is a matrix
      nans <- rep(NaN, nrow(noaa_data[[i]]))
      while (ncol(noaa_data[[i]]) < increments) {
        noaa_data[[i]] <- cbind(noaa_data[[i]], nans)
      }
    } else {   # Is a vector
      while (length(noaa_data[[i]]) < increments) {
        noaa_data[[i]] <- c(noaa_data[[i]], NaN);
      }
    }
  }
  
   #####################################
  #done with data processing- now want to take the list and make one df for upscaling to daily
  
  time = seq(from = start_date + lubridate::hours(6), to = end_date, by = "6 hour") 
  forecasts = matrix(ncol = length(noaa_data)+ 2, nrow = 0)
  colnames(forecasts) <- c(cf_var_names, "timestamp", "NOAA.member")
  
  index = matrix(ncol = length(noaa_data), nrow = length(time))
  for(i in 1:21){
    rm(index)
    index = matrix(ncol = length(noaa_data), nrow = length(time))
    for(j in 1:length(noaa_data)){
      index[,j] <- noaa_data[[j]][i,]
      colnames(index) <- c(cf_var_names) 
      index <- as.data.frame(index)
    }
    index$timestamp <- as.POSIXct(time)
    index$NOAA.member <- rep(i, times = length(time))
    forecasts <- rbind(forecasts, index)
  }
  
  forecasts <- forecasts %>% tidyr::drop_na()
  
  # convert temp from K to C
  forecasts$Temperature_height_above_ground_ens <- forecasts$Temperature_height_above_ground_ens - 273.15
  
  # change timezone and drop hours
  if(time.zone != "UTC"){
    attributes(forecasts$timestamp)$tzone <- time.zone
  }
  forecasts$timestamp <- format(forecasts$timestamp, format = "%Y-%m-%d")
  
  day.forecast <- forecasts %>% 
    group_by(NOAA.member) %>% 
    group_by(timestamp, add = TRUE) %>% 
    dplyr::summarise(max.temp = max(Temperature_height_above_ground_ens), 
                     min.temp = min(Temperature_height_above_ground_ens),
                     max.rh = max(Relative_humidity_height_above_ground_ens),
                     min.rh = min(Relative_humidity_height_above_ground_ens),
                     precip = sum(Total_precipitation_surface_6_Hour_Accumulation_ens))
  
  day.forecast <- day.forecast %>% 
    mutate(vpd = RHtoVPD(min.rh, min.temp))
  
  var.names <- c("max.temp", "min.temp", "max.rh", "min.rh", "precip", "vpd")
  cf_var_units = c("C", "C", "percent", "precent", "mm", "kPa")  
  
  n.days <- sum(day.forecast$NOAA.member==1)
  
  if (!dir.exists(outfolder)) {
    dir.create(outfolder, recursive=TRUE, showWarnings = FALSE)
  }
  

  #############################################
  # Done with data processing.  Now writing the data to the specified directory. Each ensemble member is written to its own file, for a total
  # of 21 files.  
  
  
  # Create a data frame with information about the file.  This data frame's format is an internal PEcAn standard, and is stored in the BETY database to
  # locate the data file.  The data file is stored on the local machine where the download occured.  Because NOAA GEFS is an 
  # ensemble of 21 different forecast models, each model gets its own data frame.  All of the information is the same for 
  # each file except for the file name.
  results = data.frame(
    file = "",                            #Path to the file (added in loop below).
    host = PEcAn.remote::fqdn(),          #Name of the server where the file is stored
    mimetype = "application/x-netcdf",    #Format the data is saved in
    formatname = "CF Meteorology",        #Type of data
    startdate = paste0(format(start_date, "%Y-%m-%dT%H:%M:00")),    #starting date and time, down to the second
    enddate = paste0(format(end_date, "%Y-%m-%dT%H:%M:00")),        #ending date and time, down to the second
    dbfile.name = "NOAA_GEFS_downscale",            #Source of data (ensemble number will be added later)
    stringsAsFactors = FALSE
  )
  
  results_list = list()
  
  #Each ensemble gets its own file.
  #These dimensions will be used for all 21 ncdf4 file members, so they're all declared once here.
  #The data is really one-dimensional for each file (though we include lattitude and longitude dimensions
  #to comply with the PEcAn standard).
  
  time_dim = ncdf4::ncdim_def(name="time", 
                              units = paste("Days since", format(start_date, "%Y-%m-%dT%H:%M")), 
                              1:n.days, 
                              create_dimvar = TRUE)
  # ens_dim = ncdf4::ncdim_def("Ens", "ens number", 1:21, create_dimvar = TRUE)
  lat_dim = ncdf4::ncdim_def("latitude", "degree_north", lat.in, create_dimvar = TRUE)
  lon_dim = ncdf4::ncdim_def("longitude", "degree_east", lon.in, create_dimvar = TRUE)
  
  def_list <- list()
  for(jj in seq_along(var.names)){
    def_list[[jj]] <- ncvar_def(name = var.names[jj],
                                units = cf_var_units[[jj]],
                                dim = list(time_dim, lat_dim, lon_dim),
                                missval = NaN)
  }
  
  dimensions_list = list(time_dim, lat_dim, lon_dim)
  
  nc_var_list = list()
  for (i in 1:length(var.names)) { #Each ensemble member will have data on each variable stored in their respective file.
    nc_var_list[[i]] = ncdf4::ncvar_def(var.names[i], cf_var_units[i], dimensions_list, missval=NaN)
  }
  
  day.folder = paste("NOAA_GEFS", sitename, format(start_date, "%Y-%m-%d"), sep=".")
  create.folder <- file.path(outfolder, day.folder)
  
  #Each file will go in its own folder.
  if (!dir.exists(create.folder)) {
    dir.create(create.folder, recursive=TRUE, showWarnings = FALSE)
  }

  #For each ensemble
  for (i in 1:21) { # i is the ensemble number
    #Generating a unique identifier string that characterizes a particular data set.
    identifier <- paste0(i, "_", day.folder, "_", format(end_date, "%Y-%m-%d"))
    ensemble_folder = file.path(outfolder, day.folder, identifier)
    
    data <- as.data.frame(day.forecast %>% 
                           dplyr::select(NOAA.member, var.names) %>% 
                           dplyr::filter(NOAA.member == i))
    data <- data %>% dplyr::select(var.names)  
    
    flname = paste0(ensemble_folder, ".nc")
    
    #Each ensemble member gets its own unique data frame, which is stored in results_list
    #Object references in R work differently than in other languages. When adding an item to a list, R creates a copy of it
    #for you instead of just inserting the object reference, so this works.
    results$file <- flname
    results$dbfile.name <- flname
    results_list[[i]] <- results
    
    if (!file.exists(flname) | overwrite) {
      nc_flptr = ncdf4::nc_create(flname, nc_var_list, verbose=verbose)
      
      #For each variable associated with that ensemble
      for (j in 1:length(var.names)) {
        # "j" is the variable number.  "i" is the ensemble number. 
        # Remember that each row represents an ensemble
        ncdf4::ncvar_put(nc_flptr, nc_var_list[[j]], data[,j])
      }
      
      ncdf4::nc_close(nc_flptr)  #Write to the disk/storage
    } else {
      PEcAn.logger::logger.info(paste0("The file ", flname, " already exists.  It was not overwritten."))
    }
    
  }
  file <- paste0(create.folder, "/ResultsListMeta", format(start_date, "%Y-%m-%d"), ".RData")
  save(results_list,
       file = file)
  return(results_list)
} #downscale.NOAA_GEFS

