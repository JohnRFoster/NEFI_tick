library(ecoforecastR)

source("Functions/combine_chains.R")

path <- "../FinalOut/Independent_Fits/GDDThreshold"

### WINDOW ###

subset <- "Window"

model_site <- c("GDDSwitch_GreenControl",
                "GDDSwitch_GreenControl_Cont1_",
                "GDDSwitch_HenryControl",
                "GDDSwitch_HenryControl_Cont1_",
                "GDDSwitch_TeaControl",
                "GDDSwitch_TeaControl_Cont1_",
                "Temp_GDDSwitch_GreenControl",
                "Temp_GDDSwitch_GreenControl_Cont1_",
                "Temp_GDDSwitch_HenryControl",
                "Temp_GDDSwitch_TeaControl")

for(i in 1:length(model_site)){
  load <- file.path(path,subset,model_site[i])
  save <- paste(load,".RData",sep="")
  out <- combine_chains(load, save = save)
  cat("\nCombined and saved\n",save)
}

### LowThreshOnly ###

subset <- "LowThreshOnly"

model_site <- c("GDDSwitch_GreenControl",
                "GDDSwitch_HenryControl",
                "GDDSwitch_TeaControl",
                "Temp_GDDSwitch_GreenControl",
                "Temp_GDDSwitch_HenryControl",
                "Temp_GDDSwitch_TeaControl",
                "Temp_GDDSwitch_gammaBeta33_GreenControl",
                "Temp_GDDSwitch_gammaBeta33_GreenControl_Cont1_",
                "Temp_GDDSwitch_gammaBeta33_HenryControl",
                "Temp_GDDSwitch_gammaBeta33_HenryControl_Cont1_",
                "Temp_GDDSwitch_gammaBeta33_TeaControl",
                "Temp_GDDSwitch_gammaBeta33_TeaControl_Cont1_")

for(i in 1:length(model_site)){
  load <- file.path(path,subset,model_site[i])
  save <- paste(load,".RData",sep="")
  out <- combine_chains(load, save = save)
  cat("\nCombined and saved\n",save)
}


