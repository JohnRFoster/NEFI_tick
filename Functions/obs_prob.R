

obs_prob <- function(ua, obs.temp = NULL){
  
  monitor <- names(ua)
  
  # larva obs
  if("theta.larva" %in% monitor){
    theta.larva <- ua$theta.larva
  } else {
    for(t in seq_along(obs.temp)){
      if ("beta.l.obs" %in% monitor & !("beta.l.lat" %in% monitor) & !("beta.l.vert" %in% monitor)){
        theta.larva <- 1 / (1 + ua$beta.l.obs*obs.temp[t]^2)
      } else if ("beta.l.obs" %in% monitor & "beta.l.lat" %in% monitor & !("beta.l.vert" %in% monitor)){
        theta.larva <- 1 / (1 + ua$beta.l.obs*(ua$beta.l.lat + obs.temp[t])^2)
      } else if ("beta.l.obs" %in% monitor & "beta.l.lat" %in% monitor & "beta.l.vert" %in% monitor){
        theta.larva <- 1 / (1 + ua$beta.l.vert + ua$beta.l.obs*(ua$beta.l.lat + obs.temp[t])^2)
      }
    }
  }
  
  # nymph obs
  if("theta.nymph" %in% monitor){
    theta.nymph <- ua$theta.nymph
  } else {
    for(t in seq_along(obs.temp)){
      if ("beta.n.obs" %in% monitor & !("beta.n.lat" %in% monitor) & !("beta.n.vert" %in% monitor)){
        theta.nymph <- 1 / (1 + ua$beta.n.obs*obs.temp^2)
      } else if ("beta.n.obs" %in% monitor & "beta.n.lat" %in% monitor & !("beta.n.vert" %in% monitor)){
        theta.nymph <- 1 / (1 + ua$beta.n.obs*(ua$beta.n.lat + obs.temp[t])^2)
      } else if ("beta.n.obs" %in% monitor & "beta.n.lat" %in% monitor & "beta.n.vert" %in% monitor){
        theta.nymph <- 1 / (1 + ua$beta.n.vert + ua$beta.n.obs*(ua$beta.n.lat + obs.temp[t])^2)
      }
    }
  }
  
  # adult obs
  if("theta.adult" %in% monitor){
    theta.adult <- ua$theta.adult
  } else {
    for(t in seq_along(obs.temp)){
      if ("beta.a.obs" %in% monitor & !("beta.a.lat" %in% monitor) & !("beta.a.vert" %in% monitor)){
        theta.adult <- 1 / (1 + ua$beta.a.obs*obs.temp^2)
      } else if ("beta.a.obs" %in% monitor & "beta.a.lat" %in% monitor & !("beta.a.vert" %in% monitor)){
        theta.adult <- 1 / (1 + ua$beta.a.obs*(ua$beta.a.lat + obs.temp[t])^2)
      } else if ("beta.a.obs" %in% monitor & "beta.a.lat" %in% monitor & "beta.a.vert" %in% monitor){
        theta.adult <- 1 / (1 + ua$beta.a.vert + ua$beta.a.obs*(ua$beta.a.lat + obs.temp[t])^2)
      }
    }
  }
  
  return(list(theta.larva = theta.larva,
              theta.nymph = theta.nymph,
              theta.adult = theta.adult))
}





