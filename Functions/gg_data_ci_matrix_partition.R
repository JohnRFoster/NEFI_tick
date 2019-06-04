gg_data_ci_matrix_partition <- function(CI.larvaeSurv = NULL,
                                        CI.nymphSurv = NULL,
                                        CI.adultSurv = NULL,
                                        CI.larv2nymph = NULL,
                                        CI.nymph2adult = NULL,
                                        CI.repro = NULL,
                                       life.stage, time, obs){
  
  df.out <- data.frame(time = time,
                       obs = obs)

    df.out$ci.larvaeSurv.low = CI.larvaeSurv[[life.stage]][1,] 
    df.out$ci.larvaeSurv.high = CI.larvaeSurv[[life.stage]][3,]

    df.out$ci.nymphSurv.low = CI.nymphSurv[[life.stage]][1,] 
    df.out$ci.nymphSurv.high = CI.nymphSurv[[life.stage]][3,]

    df.out$ci.adultSurv.low = CI.adultSurv[[life.stage]][1,] 
    df.out$ci.adultSurv.high = CI.adultSurv[[life.stage]][3,]

    df.out$ci.larv2nymph.low = CI.larv2nymph[[life.stage]][1,] 
    df.out$ci.larv2nymph.high = CI.larv2nymph[[life.stage]][3,]

    df.out$ci.nymph2adult.low = CI.nymph2adult[[life.stage]][1,] 
    df.out$ci.nymph2adult.high = CI.nymph2adult[[life.stage]][3,]

    df.out$ci.repro.low = CI.repro[[life.stage]][1,] 
    df.out$ci.repro.high = CI.repro[[life.stage]][3,]

  return(df.out)
}