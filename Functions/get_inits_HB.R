get_inits_HB <- function(met.proc){


  ## reduce data to just what is needed, grab inits
  if(is.null(met.proc)){
    init.dir <- "../FinalOut/A_Correct/NULL"
    load(file.path(init.dir, "Green/Combined_thinMat_NULL_GreenControl.RData"))
    init.params <- params.mat
    
    load(file.path(init.dir, "Henry/Combined_thinMat_NULL_HenryControl.RData"))
    init.params <- rbind(init.params, params.mat)
    
    load(file.path(init.dir, "Tea/Combined_thinMat_NULL_TeaControl.RData"))
    init.params <- rbind(init.params, params.mat)
    
  } else if(met.proc == "max temp"){
    
    # init.dir <- "../FinalOut/Independent_Fits/GDDThreshold/Temp_ObsProc/Obs_1_beta"
    # load(file.path(init.dir, "Green/Combined_thinMat_MaxTemp_ObsProc_beta_1_K_set_GreenControl.RData"))
    # init.params <- params.mat
    # 
    # load(file.path(init.dir, "Henry/Combined_thinMat_MaxTemp_ObsProc_beta_1_K_set_HenryControl.RData"))
    # init.params <- rbind(init.params, params.mat)
    # 
    # load(file.path(init.dir, "Tea/Combined_thinMat_MaxTemp_ObsProc_beta_1_K_set_TeaControl.RData"))
    # init.params <- rbind(init.params, params.mat)
    
  } else if(met.proc == "max rh" | met.proc == "min rh"){
    
    # init.dir <- "../FinalOut/Independent_Fits/GDDThreshold/RH_ObsProc/beta_111"
    # load(file.path(init.dir, "Green/Combined_thinMat_MaxRH_ObsProc_beta_111_K_set_GreenControl.RData"))
    # init.params <- params.mat
    # 
    # load(file.path(init.dir, "Henry/Combined_thinMat_MaxRH_ObsProc_beta_111_K_set_HenryControl.RData"))
    # init.params <- rbind(init.params, params.mat)
    # 
    # load(file.path(init.dir, "Tea/Combined_thinMat_MaxRH_ObsProc_beta_111_K_set_TeaControl.RData"))
    # init.params <- rbind(init.params, params.mat)
    
  } else if(met.proc == "vpd"){
    
    # init.dir <- "../FinalOut/Independent_Fits/GDDThreshold/VPD_ObsProc/beta_111"
    # load(file.path(init.dir, "Green/Combined_thinMat_VPD_ObsProc_beta_1_K_set_GreenControl.RData"))
    # init.params <- params.mat
    # 
    # load(file.path(init.dir, "Henry/Combined_thinMat_VPD_ObsProc_beta_1_K_set_HenryControl.RData"))
    # init.params <- rbind(init.params, params.mat)
    # 
    # load(file.path(init.dir, "Tea/Combined_thinMat_VPD_ObsProc_beta_1_K_set_TeaControl.RData"))
    # init.params <- rbind(init.params, params.mat)
    
  } else {
    stop("met.proc not set!!", call. = FALSE)
  }
  
  init.params <- apply(init.params, 2, mean)
  
  return(init.params)
}