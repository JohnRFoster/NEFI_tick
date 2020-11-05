### ================================================= ###
###      functions for getting gefs into jags         ###
### ================================================= ###

#' This function is for getting historical daily weather estimates to use as mean priors
#'
#' @param days vector of days (as day of year) to get means for
#' @example gefs_prior(c(200:210))


gefs_prior <- function(days){
  
  met <- read.csv("/projectnb/dietzelab/fosterj/Data/Cary_Met_Data_Daily.csv")
  met <- met[, c("DATE","MAX_TEMP", "MIN_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC")]
  met$DATE <- as.Date(as.character(met$DATE), format = c("%m/%d/%Y"))
  met <- met %>% 
    filter(DATE <= ymd("2017-12-31")) %>% 
    filter(DATE >= ymd("1995-01-01")) %>% 
    mutate(DOY = yday(DATE)) %>% 
    mutate(year = year(DATE))
  
  unique.years <- unique(met$year)
  
  calc_gdd <- function(base=10, max, min){
    gdd <- max(mean(max, min) - base, 0)
  }
  
  cum.gdd <- vector()
  for(i in seq_along(unique.years)){
    subset <- met %>% filter(year == unique.years[i])
    min.missing <- which(is.na(subset$MIN_TEMP))
    max.missing <- which(is.na(subset$MAX_TEMP))
    subset$MIN_TEMP[min.missing] <- mean(subset$MIN_TEMP, na.rm = TRUE)
    subset$MAX_TEMP[max.missing] <- mean(subset$MAX_TEMP, na.rm = TRUE)
    gdd <- vector()
    for(t in 1:nrow(subset)){
      gdd[t] <- calc_gdd(10, subset$MAX_TEMP[t], subset$MIN_TEMP[t]) 
    }
    cum.gdd <- c(cum.gdd, cumsum(gdd))
  }
  
  # add cumulative growing degree days to observed met
  met <- met %>% 
    as.data.frame() %>% 
    mutate(cum.gdd = cum.gdd)
  
  doy.mean <- met %>%
    group_by(DOY) %>%
    summarise(
      max.temp = mean(MAX_TEMP, na.rm = TRUE),
      min.temp = mean(MIN_TEMP, na.rm = TRUE),
      max.rh = mean(MAX_RH, na.rm = TRUE),
      min.rh = mean(MIN_RH, na.rm = TRUE),
      precip = mean(TOT_PREC, na.rm = TRUE),
      cum.gdd = mean(cum.gdd, na.rm = TRUE)
    ) %>%
    filter(DOY %in% days)
  
  return(doy.mean)
}

#' This function is for getting the mean and precision from gefs ensembles
#' for a given 16 day forecast 
#'
#' @param met.gefs list gefs ensemble, each list element is a data frame with met variables as columns
#' @param met.var met variable to extract, must match column names in met.gefs
#' 

library(LaplacesDemon)

get_gefs_mean_prec <- function(met.gefs, met.var, scale = 0){
  add.2.obs <- NULL
  dat <- t(map_dfc(met.gefs, function(x) x %>% as.data.frame() %>% select(all_of(met.var))))
  if(met.var == "cum.gdd"){
    var.check <- apply(dat, 2, var)
    if(var.check[1] == 0){
      dat <- dat[,-1]
      add.2.obs <- dat[1,1] # only change when assessing cum.gdd and var = 0
    } 
  } else {
    dat <- dat - scale
  }
  
  mu <- apply(dat, 2, mean) # mean vector
  prec.mat <- solve(cov(dat), tol = 1e-30) # precision matrix
  return(list(mu = mu, prec = prec.mat, add.2.obs = add.2.obs))
}

#' This function reads the files for a given gefs forecast and creates a list
#' each element in the list is a data frame for one ensemble forecast
#'
#' @param dir file path to directory ensembles are stored 
#' @param end.cum.gdd cumulative growing degree days at the end of met observations
#' 
#' # directories

get_gefs_rnoaa <- function(date, end.cum.gdd){
  
  dir.gefs.rnoaa <- "/projectnb/dietzelab/fosterj/Data/GEFSrnoaa/Cary/" # from rnoaa downloads
  day.dir <- paste0("NOAA_GEFS.CaryInstitute.", date) # file name with date 
  ens.store <- paste0(dir.gefs.rnoaa, day.dir) # dir with individual ensembles
  
  ens.files <- list.files(ens.store)
  ens.files <- ens.files[grepl(".nc", ens.files)] # want only .nc files
  
  # get dimensions
  gefs.path <- file.path(ens.store, ens.files[1]) # full path to .nc 
  gefs.ens <- nc_open(gefs.path) 
  n.var <- length(gefs.ens$var)
  var.names <- names(gefs.ens$var)
  n.days <- length(ncvar_get(gefs.ens, var.names[1]))
  
  met.gefs <- list()
  for(ens in 1:21){ # get data for each ensemble member
    gefs.path <- file.path(ens.store, ens.files[ens]) # full path to .nc 
    gefs.ens <- nc_open(gefs.path)  
    gefs.data <- matrix(NA, n.days, n.var)
    for(v in seq_along(var.names)){
      gefs.data[,v] <- ncvar_get(gefs.ens, var.names[v])  
    }
    colnames(gefs.data) <- var.names
    
    gdd <- rep(NA, nrow(gefs.data))
    for(i in 1:nrow(gefs.data)){
      gdd[i] <- max(mean(gefs.data[i, "max.temp"], gefs.data[i, "min.temp"]) - 10, 0)
    }
    
    gefs.cum.gdd <- cumsum(c(end.cum.gdd, gdd))
    gefs.data <- gefs.data %>% 
      as.data.frame() %>% 
      mutate(cum.gdd = gefs.cum.gdd[-1])
    
    met.gefs[[ens]] <- gefs.data
  }
  
  # the returned object is a list 
  # each element represents an ensemble member
  # each element is a data frame, days in rows, vars in columns
  return(met.gefs)
}

get_gefs_noaaGEFSpoint <- function(method, date, end.cum.gdd){
  
  # directories by method
  if(method == "grid"){
    dir.gefs <- "/projectnb/dietzelab/fosterj/Data/GEFSgrid/NOAAGEFS_6hr/CARY/" # from noaaGEFSgrid
    day.dir <- as.character(date) # dir name with date, grid
  } else if(method == "point"){
    dir.gefs <- "/projectnb/dietzelab/fosterj/Data/GEFSpoint/NOAAGEFS_6hr/CARY/" # from noaaGEFSpoint
    day.dir <- as.character(gsub("-", "", date)) # dir name with date, point
  }
  
  ens.store <- paste0(dir.gefs, day.dir, "/00") # dir with individual ensembles
  ens.files <- list.files(ens.store)
  ens.files <- ens.files[grepl(".nc", ens.files)] # want only .nc files
  
  met.gefs <- list()
  for(ens in seq_along(ens.files)){
    gefs.path <- file.path(ens.store, ens.files[ens]) # full path to .nc 
    gefs.ens <- nc_open(gefs.path)
    
    time.step <- ncvar_get(gefs.ens, "time") # hours since forecast date
    time <- date + (time.step / 24) # forecast date
    tz(time) <- "UTC" # forecasts are UTC
    time <- with_tz(time, "America/New_York") # change to local time, keeps instance
    
    # get variables
    temp <- ncvar_get(gefs.ens, "air_temperature") # units = K
    rh <- ncvar_get(gefs.ens, "relative_humidity") # units = percent
    # precip <- ncvar_get(gefs.ens, "precipitation_flux") # units = kg m-2 s-1

    met.df <- data.frame(time = time,
                         temp = temp - 273.15, # units = C
                         rh = rh) # units = percent
                         # precip = precip)
    
    gefs.daily <- met.df %>% 
      mutate(day = floor_date(time, "day")) %>%
      group_by(day) %>% 
      summarise(max.temp = max(temp),
                min.temp = min(temp),
                max.rh = max(rh),
                min.rh = min(rh)) %>% 
      filter(day >= date)
    
    gdd <- rep(NA, nrow(gefs.daily))
    for(i in 1:nrow(gefs.daily)){
      gdd[i] <- max(mean(pull(gefs.daily[i, "max.temp"]),
                         pull(gefs.daily[i, "min.temp"])) - 10, 0)
    }
    
    gefs.cum.gdd <- cumsum(c(end.cum.gdd, gdd))
    gefs.daily <- gefs.daily %>% 
      mutate(cum.gdd = gefs.cum.gdd[-1]) %>% 
      select(-day)
    
    # some ensembles from grid method only have 16 days, 
    # so lets get rid of those and keep the 35 day forecasts
    if(method == "grid"){
      if(nrow(gefs.daily) < 17) gefs.daily <- NULL
    }
    
    met.gefs[[ens]] <- gefs.daily
  }
  
  # the returned object is a list 
  # each element represents an ensemble member
  # each element is a data frame, days in rows, vars in columns
  met.gefs <- purrr::compact(met.gefs)
  return(met.gefs)
}

get_gefs_archive <- function(date, end.cum.gdd){
  dir.gefs.arch <- "/projectnb/dietzelab/fosterj/Data/GEFSarchive/" # archived GEFS downloads
  pattern <- as.character(gsub("-", "", date)) # beginning of csv
  
  ens.files <- list.files(dir.gefs.arch)
  ens.files <- ens.files[grepl(".csv", ens.files)] # want only .csv files

  csv.index <- which(str_detect(ens.files, pattern)) # index of file we want

  gefs.csv <- read.csv(paste0(dir.gefs.arch, ens.files[csv.index]))  # read file
  gefs.csv$forecast.date <- ymd_hms(gefs.csv$forecast.date)
  tz(gefs.csv$forecast.date) <- "GMT" # forecasts are GMT, set in process script
  gefs.csv$forecast.date <- with_tz(gefs.csv$forecast.date, "America/New_York") # change to local time, keeps instance
  
  # number of ensembles
  n.ens <- gefs.csv %>% 
    pull(ensembles) %>% 
    unique() %>% 
    length()
  
  # subset and filter 
  gefs.daily <- gefs.csv %>% 
    select(all_of(c("forecast.date", "ensembles", "tmp2m", "rh2m"))) %>% 
    mutate(day = floor_date(forecast.date, "day")) %>%
    mutate(tmp2m = tmp2m - 273.15) %>%
    group_by(ensembles, day) %>% 
    summarise(max.temp = max(tmp2m),
              min.temp = min(tmp2m),
              max.rh = max(rh2m),
              min.rh = min(rh2m)) %>% 
    filter(day >= date)
  
  met.gefs <- list()
  for(ens in 1:n.ens){
    gefs.subset <- gefs.daily %>% 
      filter(ensembles == ens)
    
    gdd <- rep(NA, nrow(gefs.subset))
    for(i in 1:nrow(gefs.subset)){
      gdd[i] <- max(mean(pull(gefs.subset[i, "max.temp"]), 
                         pull(gefs.subset[i, "min.temp"])) - 10, 0)
    }
    
    gefs.cum.gdd <- cumsum(c(end.cum.gdd, gdd))
    gefs.subset <- gefs.subset %>% 
      mutate(cum.gdd = gefs.cum.gdd[-1]) %>% 
      select(-day)
    
    
    
    met.gefs[[ens]] <- gefs.subset
  }
  
  
  
  # the returned object is a list 
  # each element represents an ensemble member
  # each element is a data frame, days in rows, vars in columns
  return(met.gefs)
  
}


