library(rNOMADS)
library(tidyverse)
library(lubridate)


Sys.setenv(PATH = "/projectnb/dietzelab/fosterj/grib2/wgrib2/")
nens <- 21
top.dir <- "/projectnb/dietzelab/fosterj/"
working_directory <- paste0(top.dir, "Data/GEFSarchive/")

day_vec <- as.character(20200622)

# read array number and subset day_vec
# array.num <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# day_list <- day_vec[array.num]
day_list <- day_vec

cat("Extracting", day_list, "\n\n")

#ISSUE WITH 20181223 gens-a_3_20181223_1200_054_14.grb2

ensemble <- c('00','01','02','03','04','05','06','07','08','09','10',
              '11','12','13','14','15','16','17','18','19','20')

start_time <- '1200'

hours <- c('000','006','012','018','024','030','036','042','048','054','060','066','072','078',
           '084','090','096','102','108','114','120','126','132','138','144','150','156','162',
           '168','174','180','186','192','198','204','210','216','222','228','234','240','246',
           '252','258','264','270','276','282','288','294','300','306','312','318','324','330',
           '336','342','348','354','360','366','372','378','384')

lake_lat_n <- 41.7851
lake_lon_w <- 73.7338
lon.dom <- seq(0, 360, by = 1) #domain of longitudes in model
lat.dom <- seq(-90, 90, by = 1) #domain of latitudes in model
lon <- which.min(abs(lon.dom  - (360 - lake_lon_w))) - 1 #NOMADS indexes start at 0
lat <- which.min(abs(lat.dom - lake_lat_n)) - 1 #NOMADS indexes start at 0 

for(d in 1:length(day_list)){
  
  AirTemp <- array(9.999e+20,dim=c(length(hours),nens))
  LongWave <- array(9.999e+20,dim=c(length(hours),nens))
  ShortWave <- array(9.999e+20,dim=c(length(hours),nens))
  RelHum <- array(9.999e+20,dim=c(length(hours),nens))
  WindSpeed <- array(9.999e+20,dim=c(length(hours),nens))
  ugrd10m <- array(9.999e+20,dim=c(length(hours),nens))
  vgrd10m <- array(9.999e+20,dim=c(length(hours),nens))
  Rain <- array(9.999e+20,dim=c(length(hours),nens))
  Snow <- array(9.999e+20,dim=c(length(hours),nens))
  tcdcclm <- array(9.999e+20,dim=c(length(hours),nens))
  pressfc <- array(9.999e+20,dim=c(length(hours),nens))
  spfh2m <- array(9.999e+20,dim=c(length(hours),nens))
  ensembles <- array(9.999e+20,dim=c(length(hours),nens))
  forecast.date <- array(9.999e+20,dim=c(length(hours),nens))

  day_directory <- paste0(working_directory,day_list[d],"12/")
  day_directory <- paste0(working_directory,day_list[d], "/")
  
  for(ens in 1:nens){

    
    ensemble_name <- ensemble[ens]
    if(day_list[d] == '20181221' & ens == 15){
      ensemble_name <- ensemble[14]
    }
    if(day_list[d] == '20181221' & ens == 9){
      ensemble_name <- ensemble[8]
    }
    if(day_list[d] == '20181221' & ens == 16){
      ensemble_name <- ensemble[17]
    }
    if(day_list[d] == '20181222' & ens == 12){
      ensemble_name <- ensemble[13]
    }
    
    if(day_list[d] == '20181223' & ens == 15){
      ensemble_name <- ensemble[14]
    }
    
    if(day_list[d] == '20181226' & ens == 4){
      ensemble_name <- ensemble[5]
    }
    
    if(day_list[d] == '20181227' & ens == 12){
      ensemble_name <- ensemble[11]
    }
    
    if(day_list[d] == '20181228' & ens == 7){
      ensemble_name <- ensemble[5]
    }
    
    ens_base_a <- paste(day_directory,"gens_3_",day_list[d],'12_',ensemble_name,'.g2/gens-a_3_',day_list[d],sep='')
    ens_base_b <- paste(day_directory,"gens_3_",day_list[d],'12_',ensemble_name,'.g2/gens-b_3_',day_list[d],sep='')
    
    tar_file <- paste0(day_directory,"gens_3_",day_list[d],'12_',ensemble_name,".g2.tar")
    tar_location <- paste0(day_directory,"gens_3_",day_list[d],'12_',ensemble_name,".g2")
    untar(tar_file,exdir = tar_location)
    for(hr in 1:length(hours)){
      print(c(day_list[d], ens, hr))
      
      
      gribfile_a <- paste(ens_base_a,'_',start_time,'_',hours[hr],'_',ensemble_name,'.grb2',sep='')
      gribfile_b <- paste(ens_base_b,'_',start_time,'_',hours[hr],'_',ensemble_name,'.grb2',sep='')
      if(hr == 1){
        grib_a <- ReadGrib(gribfile_a, levels = c("2 m above ground", 
                                                  "10 m above ground",
                                                  "surface"), 
                           variables = c("TMP","RH", "VGRD", "UGRD", "PRES"))
        grib_b <- ReadGrib(gribfile_b, levels = c("2 m above ground"), 
                           variables = c("SPFH"))  
      }else{
        
        grib_a <- ReadGrib(gribfile_a, levels = c("2 m above ground", 
                                                  "10 m above ground",
                                                  "surface",
                                                  "entire atmosphere"), 
                           variables = c("TMP","RH", "VGRD", "UGRD", "PRES","TCDC","DSWRF","DLWRF"))
        grib_b <- ReadGrib(gribfile_b, levels = c("2 m above ground","surface"), 
                           variables = c("SPFH","PRATE"))
      }
      
      
      grib_a_grid <- ModelGrid(grib_a, c(1.0, 1.0))
      grib_b_grid <- ModelGrid(grib_b, c(1.0, 1.0))
      
      forecast.date[hr,ens] <- as.numeric(grib_a_grid$fcst.date)
      ensembles[hr,ens] <- ens
      
      index1 <- which(grib_a_grid$levels == "2 m above ground")
      index2 <- which(grib_a_grid$variables == "TMP")
      AirTemp[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
      
      index1 <- which(grib_a_grid$levels == "2 m above ground")
      index2 <- which(grib_a_grid$variables == "RH")
      RelHum[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
      
      index1 <- which(grib_a_grid$levels == "10 m above ground")
      index2 <- which(grib_a_grid$variables == "VGRD")
      vgrd10m[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
      
      index1 <- which(grib_a_grid$levels == "10 m above ground")
      index2 <- which(grib_a_grid$variables == "UGRD")
      ugrd10m[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
      
      index1 <- which(grib_a_grid$levels == "surface")
      index2 <- which(grib_a_grid$variables == "PRES")
      pressfc[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
      
      
      index1 <- which(grib_b_grid$levels == "2 m above ground")
      index2 <- which(grib_b_grid$variables == "SPFH")
      spfh2m[hr,ens] <- grib_b_grid$z[index1, index2 , lon,lat]
      
      if(hr > 1){
        
        index1 <- which(grib_a_grid$levels == "entire atmosphere")
        index2 <- which(grib_a_grid$variables == "TCDC")
        tcdcclm[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
        
        index1 <- which(grib_a_grid$levels == "surface")
        index2 <- which(grib_a_grid$variables == "DSWRF")
        ShortWave[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
        
        index1 <- which(grib_a_grid$levels == "surface")
        index2 <- which(grib_a_grid$variables == "DLWRF")
        LongWave[hr,ens] <- grib_a_grid$z[index1, index2 , lon,lat]
        
        index1 <- which(grib_b_grid$levels == "surface")
        index2 <- which(grib_b_grid$variables == "PRATE")
        Rain[hr,ens] <- grib_b_grid$z[index1, index2 , lon,lat]
      }
    }
    unlink(tar_location, recursive = TRUE)
  }
  
  forecast.date <- as.numeric(forecast.date)
  
  forecast.time <- as_datetime(forecast.date, origin = lubridate::origin,tz = 'GMT')
  
  forecast_noaa <- data.frame(forecast.date = forecast.time, 
                              ensembles = c(ensembles), 
                              tmp2m = c(AirTemp), 
                              dlwrfsfc= c(LongWave), 
                              dswrfsfc = c(ShortWave), 
                              pratesfc = c(Rain), 
                              rh2m = c(RelHum), 
                              vgrd10m = c(vgrd10m), 
                              ugrd10m = c(ugrd10m), 
                              spfh2m = c(spfh2m), 
                              pressfc = c(pressfc),
                              tcdcclm = c(tcdcclm))
  
  file_name <- paste0(working_directory,day_list[d],"gep_all_12z.csv")
  write.csv(forecast_noaa,file_name,row.names = FALSE)
}