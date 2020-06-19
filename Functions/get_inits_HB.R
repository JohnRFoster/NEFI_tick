get_inits_HB <- function(init.dir = NULL, init.model.name = NULL, met.proc = NULL){


  ## reduce data to just what is needed, grab inits
  if(is.null(met.proc)){
    init.dir <- "../FinalOut/A_Correct/NULL"
    load(file.path(init.dir, "Green/Combined_thinMat_NULL_GreenControl.RData"))
    init.params <- params.mat
    
    load(file.path(init.dir, "Henry/Combined_thinMat_NULL_HenryControl.RData"))
    init.params <- rbind(init.params, params.mat)
    
    load(file.path(init.dir, "Tea/Combined_thinMat_NULL_TeaControl.RData"))
    init.params <- rbind(init.params, params.mat)
    
    init.params <- apply(init.params, 2, mean)
    return(init.params)
    
  } 
  
  if(!is.null(init.dir) & !is.null(init.model.name)){
    load(paste0(init.dir, "Green/", init.model.name, "_GreenControl.RData"))
    init.params <- params.mat
    
    load(paste0(init.dir, "Henry/", init.model.name, "_HenryControl.RData"))
    init.params <- rbind(init.params, params.mat)
    
    load(paste0(init.dir, "Tea/", init.model.name, "_TeaControl.RData"))
    init.params <- rbind(init.params, params.mat)
    
    init.params <- apply(init.params, 2, mean)
    return(init.params)
  }
  
  
}