get_last_day <- function(){

  sites <- c("Green Control","Henry Control","Tea Control")
  raw.dat <- read.csv("tick_cleaned")   # read in data
  raw.dat$DATE <- as.character(raw.dat$DATE) # convert to date
  
  # find last day in timeseries for each site
  n_obs <- last.day <- vector()
  for(i in 1:length(sites)){
    subset <- subset(raw.dat, Grid == sites[i])
    n_obs[i] <- nrow(subset)
    last.day[i] <- as.character(subset$DATE[n_obs[i]])
  }
  
return(data.frame(n.obs = n_obs,
                  last.day = last.day))
}

