gg_data_ci_partition <- function(CI.ic = NULL,
                                 CI.param = NULL,
                                 CI.process = NULL,
                                 CI.ic.param.proc = NULL,
                                 CI.ic.param = NULL,
                                 CI.all = NULL,
                                 CI.random = NULL,
                                 median = NULL,
                                 life.stage, time, obs){
  
  df.out <- data.frame(time = time,
                       obs = obs,
                       median = median)
  if(!is.null(CI.ic)){
    df.out$ci.ic.low = CI.ic[[life.stage]][1,] 
    df.out$ci.ic.high = CI.ic[[life.stage]][3,]
  }
  if(!is.null(CI.param)){
    df.out$ci.param.low = CI.param[[life.stage]][1,] 
    df.out$ci.param.high = CI.param[[life.stage]][3,]
  }
  if(!is.null(CI.process)){
    df.out$ci.process.low = CI.process[[life.stage]][1,] 
    df.out$ci.process.high = CI.process[[life.stage]][3,]
  }
  if(!is.null(CI.all)){
    df.out$ci.all.low = CI.all[[life.stage]][1,] 
    df.out$ci.all.high = CI.all[[life.stage]][3,]
  }
  if(!is.null(CI.random)){
    df.out$ci.random.low = CI.random[[life.stage]][1,] 
    df.out$ci.random.high = CI.random[[life.stage]][3,]
  }
  if(!is.null(CI.ic.param.proc)){
    df.out$ci.ic.param.proc.low = CI.ic.param.proc[[life.stage]][1,] 
    df.out$ci.ic.param.proc.high = CI.ic.param.proc[[life.stage]][3,]
  }
  if(!is.null(CI.ic.param)){
    df.out$ci.ic.param.low = CI.ic.param[[life.stage]][1,] 
    df.out$ci.ic.param.high = CI.ic.param[[life.stage]][3,]
  }

  return(df.out)
}