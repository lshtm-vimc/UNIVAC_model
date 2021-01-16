#Combine univac PSA Input Paramters output file for each country into one csv file
rm(list = ls())

pathogen_list = c("Rota", "Sp", "Hib")
pathogen = pathogen_list[2]
wd = paste0("~/UNIVAC_model/Sandbox/Output/VIMC_", pathogen, "_Results/PSA")
setwd(wd)
require(data.table)
require(tidyverse)


files_to_df <- function(path, pattern){ 
  
  # pattern <- "data"
  myfiles <- list.files(path = path, recursive = T, pattern = pattern, full.names = T) 
  
  df_list <- lapply(myfiles, read.csv, header = TRUE) 
  
  # Name each dataframe with the run and filename
  #names(df_list) <- str_extract(myfiles, "run\\d+")
  
  # Create combined dataframe  
  df <- df_list %>%
    bind_rows(.id = 'run')
  
  #sort by run, country and age
  df = data.table(df)
  df = setorderv(df, c("country"))
  
  # Assign dataframe to the name of the pattern  
  assign(pattern, df)
  
  # Return the dataframe  
  return(df)
  
}


ddmm <- list.dirs(path = ".", full.names = F, recursive = F)
# if(pathogen == "Hib"){
#   novax_list = ddmm[6:9]
# }else{
#   novax_list = ddmm
# }

n = 0
for(i in novax_list){
  n = n + 1
  path = i
  df <- files_to_df(path = path, pattern = "InputParams")
  assign(paste0("InputParams_df", n), df)
  
}

rm(df)
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
InputParams_all = do.call(rbind, dfs)
InputParams_all = InputParams_all[,-1]
setcolorder(InputParams_all, c("run.id", "country", colnames(InputParams_all)[3:ncol(InputParams_all)]))
colnames(InputParams_all)[1] <- "run_id"
InputParams_all = unique(InputParams_all)

fwrite(InputParams_all, file = paste0("VIMC_", pathogen, "_stochastic_params.csv"), sep = ",", col.names = T)
 
