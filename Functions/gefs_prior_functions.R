### ================================================= ###
###      functions for getting gefs into jags         ###
### ================================================= ###

#' This function is for getting historical daily weather estimates to use as mean priors
#'
#' @param days vector of days (as day of year) to get means for
#' @example gefs_prior(c(200:210))


gefs_prior <- function(days){
  
  met <- read.csv("../Cary_Met_Data_Daily.csv")
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
  prec.mat <- solve(cov(dat)) # precision matrix
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

get_gefs_point <- function(date, end.cum.gdd){
  dir.gefs.point <- "/projectnb/dietzelab/fosterj/Data/GEFSpoint/NOAAGEFS_6hr/CARY/" # from noaaGEFSpoint
  # will have to adjust from utc to eastern
}

get_gefs_grid <- function(date, end.cum.gdd){
  dir.gefs.grid <- "/projectnb/dietzelab/fosterj/Data/GEFSgrid/NOAAGEFS_6hr/CARY/" # from noaaGEFSgrid
  # will have to adjust from utc to eastern
}

get_gefs_archive <- function(date, end.cum.gdd){
  dir.gefs.arch <- "/projectnb/dietzelab/fosterj/Data/GEFSarchive/" # archived GEFS downloads
  # will have to adjust from utc to eastern
}


