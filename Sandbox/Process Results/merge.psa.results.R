#Combine univac output files for each country into one csv file
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
  names(df_list) <- str_extract(myfiles, "run\\d+")
  
  # Create combined dataframe  
  df <- df_list %>%
    bind_rows(.id = 'run')
  
  #sort by run, country and age
  df = data.table(df)
  df = setorderv(df, c("country", "age", "run_id"))
  
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
  path = paste0(i, "/no_vax")
  df <- files_to_df(path = path, pattern = "NoVax")
  assign(paste0("NoVax_df", n), df)
  
}

rm(df)
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
NoVax_all = do.call(rbind, dfs)
NoVax_all_order = setorderv(NoVax_all, c("country", "run_id"))
NoVax_all_order = NoVax_all_order[,-1]
NoVax_all_order = unique(NoVax_all_order)
if(pathogen == "Sp"){
  NoVax_all_order$disease = "PCV"
}
fwrite(NoVax_all_order, file = paste0("VIMC_", pathogen, "_PSA_NoVax.csv"), sep = ",", col.names = T)

rm(list=setdiff(ls(pattern="NoVax"), lsf.str()))

n = 0
for(i in novax_list){
  n = n + 1
  path = paste0(i, "/default_routine")
  df <- files_to_df(path = path, pattern = "default")
  assign(paste0("default_df", n), df)
  
}

rm(df)
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
df_all = do.call(rbind, dfs)
df_all_order = setorderv(df_all, c("country", "run_id"))
df_all_order = df_all_order[,-1]
df_all_order = unique(df_all_order)
if(pathogen == "Sp"){
  df_all_order$disease = "PCV"
}
fwrite(df_all_order, file = paste0("VIMC_", pathogen, "_PSA_Default.csv"), sep = ",", col.names = T)

rm(list=setdiff(ls(pattern="df"), lsf.str()))

n = 0
for(i in ddmm){
  n = n + 1
  path = paste0(i, "/best_routine")
  df <- files_to_df(path = path, pattern = "best")
  assign(paste0("best_df", n), df)
  
}

rm(df)
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
df_all = do.call(rbind, dfs)
df_all_order = setorderv(df_all, c("country", "run_id"))
df_all_order = df_all_order[,-1]
df_unique = unique(df_all_order)
if(pathogen == "Sp"){
  df_unique$disease = "PCV"
}
fwrite(df_unique, file = paste0("VIMC_", pathogen, "_PSA_Best.csv"), sep = ",", col.names = T)


