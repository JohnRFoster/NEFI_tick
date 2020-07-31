##' @title Download Data from a NOAA met station
##' 
##'
##' @param outfolder Directory where results should be written
##' @param sitename The unique ID given to each site. This is used as part of the file name.
##' @param overwrite logical. Download a fresh version even if a local file with the same name already exists?
##' @param station.id The station to download from
##' @export
##' 
##' 
##'

library(stringr)


NOAA_observed <- function(outfolder, sitename, station.id, 
                          dataset.id = "GHCND", overwrite = FALSE, day.2.grab = NULL){
  
  #Each file will go in its own folder.
  if (!dir.exists(outfolder)) {
    dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  }
  
  noaa.api.key <- "DsPttxqoSKEPGpkLKdlNOYLhyTvLUfDC"
  
  # check available data
  data.avail <- ncdc_datatypes(stationid = station.id,
                               datasetid = dataset.id,
                               token = noaa.api.key)$data
  
  if(is.null(day.2.grab)){
    fls <- list.files(outfolder) # days already downloaded
    last <- fls[length(fls)] # last file downloaded
    last.date <- as.Date(str_extract(last, "\\d{4}-\\d{2}-\\d{2}"), "%Y-%m-%d") # last day
    
    # most recent observed data minus one day to make sure there is data 
    # sometimes most recent day is empty
    grab.date <- ymd(data.avail$maxdate) - 1 
    grab.date <- as.Date(min(grab.date), "%Y-%m-%d")  # make sure we use same day for all data
    
    if(grab.date == (last.date + 1)){ # days are in sequence
      cat("Last day downloaded", format(last.date, "%Y-%m-%d"), "\n")
      cat("Downloading", format(grab.date, "%Y-%m-%d"), "\n")
    } else { # if they are not in sequence
      grab.date <- seq.Date(last.date+1, grab.date, by = 1)
      cat("There are missing days. Downloading", format(grab.date, "%Y-%m-%d"))
    }
  } else {
    grab.date <- ymd(day.2.grab)
  }
  
  for(i in seq_along(grab.date)){
    out <- ncdc(
      datasetid = dataset.id,
      datatypeid = data.avail$id,
      stationid = station.id,
      startdate = format(grab.date[i], "%Y-%m-%d"),
      enddate = format(grab.date[i], "%Y-%m-%d"),
      add_units = TRUE,
      token = noaa.api.key
    )$data
    
    # units are in tenths, so need to being back to whole units
    out$value <- out$value / 10 
    units <- c("mm", "celcius", "celcius")
    var.names <- c("PRCP", "TMAX", "TMIN")
    
    observed.prcp <- out %>% 
      filter(datatype == "PRCP") %>% 
      select(value) %>% 
      as.numeric()
    observed.tmax <- out %>% 
      filter(datatype == "TMAX") %>% 
      select(value) %>% 
      as.numeric()
    observed.tmin <- out %>% 
      filter(datatype == "TMIN") %>% 
      select(value) %>% 
      as.numeric()
    
    prcp_dim <- ncdf4::ncdim_def(name = "precipitation",
                                 units = "mm",
                                 vals = 1,
                                 create_dimvar = TRUE)
    tmax_dim <- ncdf4::ncdim_def(name = "maxTemperature",
                                 units = "celcius",
                                 vals = 1,
                                 create_dimvar = TRUE)
    tmin_dim <- ncdf4::ncdim_def(name = "minTemperature",
                                 units = "celcius",
                                 vals = 1,
                                 create_dimvar = TRUE)
    
    def_list <- list()
    for(jj in seq_along(var.names)){
      def_list[[jj]] <- ncvar_def(name = var.names[jj],
                                  units = units[jj],
                                  dim = list(prcp_dim, 
                                             tmax_dim,
                                             tmin_dim),
                                  missval = NaN)
    }
    
    day.file <- paste("NOAA_Observed", sitename, format(grab.date[i], "%Y-%m-%d"), sep = ".")
    create.file <- paste0(file.path(outfolder, day.file), ".nc")
    
    if (!file.exists(create.file) | overwrite) {
      ncout <- nc_create(create.file, def_list)
      
      ncvar_put(ncout, def_list[[1]], observed.prcp)
      ncvar_put(ncout, def_list[[2]], observed.tmax)
      ncvar_put(ncout, def_list[[3]], observed.tmin)
      
      # Global file metadata
      ncatt_put(ncout,0,"ObservationDate", as.character(format(grab.date, "%Y-%m-%d")), 
                prec =  "text")
      ncatt_put(ncout,0,"StationID",as.character(station.id), 
                prec =  "text")
      nc_close(ncout)
      
    } else {
      cat("The file ", create.file, "already exists.  It was not overwritten. \n")
    }  
  }
}