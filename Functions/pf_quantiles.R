# this function takes the output from the particle filter and
# builds the quantiles for each day from the 16-day hindcasat



pf_quantiles <- function(hindcast_out, obs.index, Nmc){
  # start counter
  obs.day <- 0 
  
  # storage
  hind.out <- matrix(NA, 5, dim(hindcast_out)[3]-1)
  ens.observation.days <- list() 
  build.index <- c(1, obs.index+1, dim(hindcast_out)[3])
  
  for(t in 1:(length(build.index)-1)){
    # days sequence
    start <- build.index[t]
    end <- build.index[t+1]-1
    
    # loop over days
    for(day in start:end){
      if(day == start){ # first day
        hind.ens <- hindcast_out[1,,day]
      } else if(day < start+16){ # before 16 days out
        day.seq <- (day-1):start
        ens.seq <- 1:length(day.seq)
        hind.ens <- matrix(NA, length(ens.seq), Nmc)
        for(j in ens.seq){
          hind.ens[j,] <- hindcast_out[j,,day.seq[j]]
        }
      } else { # after 16 days out
        day.seq <- (day-1):(day-16)
        ens.seq <- 1:16
        hind.ens <- matrix(NA, 16, Nmc)
        for(j in ens.seq){
          hind.ens[j,] <- hindcast_out[j,,day.seq[j]]
        }
      }
      hind.out[,day] <- quantile(hind.ens, c(0.025, 0.25, 0.5, 0.75, 0.975))
      
      # counter adds 1 when at observation day
      # store PF forecast for observation day
      if(day %in% obs.index){
        obs.day <- obs.day + 1
        ens.observation.days[[obs.day]] <- hind.ens
      }
    }
  }
  return(list(hind.out.quantile = hind.out,
              ens.observation.days = ens.observation.days))
}