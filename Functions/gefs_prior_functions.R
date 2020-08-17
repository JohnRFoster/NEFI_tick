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

get_gefs_mean_prec <- function(met.gefs, met.var){
  dat <- t(map_dfc(met.gefs, function(x) x %>% as.data.frame() %>% select(all_of(met.var))))
  mu <- apply(dat, 2, mean) # mean vector
  prec.mat <- solve(cov(dat), tol = 1e-30) # precision matrix
  return(list(mu = mu, prec = prec.mat))
}
