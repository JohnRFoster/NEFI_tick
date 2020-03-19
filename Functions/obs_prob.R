

obs_prob <- function(ua, N_est, obs.temp = NULL){
  
  monitor <- names(ua)
  obs.driver <- obs.temp
  
  # larva obs
  if("theta.larva" %in% monitor){
    theta.larva <- function(){ua$theta.larva}
  } else if ("beta.l.obs" %in% monitor & !("beta.l.lat" %in% monitor) & !("beta.l.vert" %in% monitor)){
    theta.larva <- function(){1 / (1 + ua$beta.l.obs*obs.temp^2)}
  } else if ("beta.l.obs" %in% monitor & "beta.l.lat" %in% monitor & !("beta.l.vert" %in% monitor)){
    theta.larva <- function(){1 / (1 + ua$beta.l.obs*(ua$beta.l.lat + obs.temp)^2)}
  } else if ("beta.l.obs" %in% monitor & "beta.l.lat" %in% monitor & "beta.l.vert" %in% monitor){
    theta.larva <- function(){1 / (1 + ua$beta.l.vert + ua$beta.l.obs*(ua$beta.l.lat + obs.temp)^2)}
  }
  
  # nymph obs
  if("theta.nymph" %in% monitor){
    theta.nymph <- function(){ua$theta.nymph}
  } else if ("beta.n.obs" %in% monitor & !("beta.n.lat" %in% monitor) & !("beta.n.vert" %in% monitor)){
    theta.nymph <- function(){1 / (1 + ua$beta.n.obs*obs.temp^2)}
  } else if ("beta.n.obs" %in% monitor & "beta.n.lat" %in% monitor & !("beta.n.vert" %in% monitor)){
    theta.nymph <- function(){1 / (1 + ua$beta.n.obs*(ua$beta.n.lat + obs.temp)^2)}
  } else if ("beta.n.obs" %in% monitor & "beta.n.lat" %in% monitor & "beta.n.vert" %in% monitor){
    theta.nymph <- function(){1 / (1 + ua$beta.n.vert + ua$beta.n.obs*(ua$beta.n.lat + obs.temp)^2)}
  }
  
  # adult obs
  if("theta.adult" %in% monitor){
    theta.adult <- function(){ua$theta.adult}
  } else if ("beta.a.obs" %in% monitor & !("beta.a.lat" %in% monitor) & !("beta.a.vert" %in% monitor)){
    theta.adult <- function(){1 / (1 + ua$beta.a.obs*obs.temp^2)}
  } else if ("beta.a.obs" %in% monitor & "beta.a.lat" %in% monitor & !("beta.a.vert" %in% monitor)){
    theta.adult <- function(){1 / (1 + ua$beta.a.obs*(ua$beta.a.lat + obs.temp)^2)}
  } else if ("beta.a.obs" %in% monitor & "beta.a.lat" %in% monitor & "beta.a.vert" %in% monitor){
    theta.adult <- function(){1 / (1 + ua$beta.a.vert + ua$beta.a.obs*(ua$beta.a.lat + obs.temp)^2)}
  }
  
  larva.prob <- nymph.prob <- adult.prob <- matrix(NA, length(ua$deviance), N_est)
  for(t in 1:N_est){
    obs.temp <- obs.driver[t]
    larva.prob[,t] <- theta.larva()
    nymph.prob[,t] <- theta.nymph()
    adult.prob[,t] <- theta.adult()
  }
  
  return(list(theta.larva = larva.prob,
              theta.nymph = nymph.prob,
              theta.adult = adult.prob))
}





