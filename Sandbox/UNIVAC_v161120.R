
## Static Model R version of the UNIVAC Excel Model
## Authors: Hira Tanvir, Andrew Clark 
## Date: 16 November 2020 -- latest version 

rm(list=ls())
source('univac_functions_v161120.R')

library(tictoc)
require(pbapply)
library(data.table)
library(extraDistr)
require(foreach)
require(doParallel)
library(mypackage)
require("parallel")

#country-specific inputs
mid.pop            = fread("Input/Pop_Mid.csv", header = TRUE)
low.pop            = fread("Input/Pop_Low.csv", header = TRUE)
high.pop           = fread("Input/Pop_High.csv", header = TRUE)
u5mr_dat           = fread("Input/vimc_u5mr.csv",header = T)  
rates_all          = fread("Input/vimc_rates.csv", header = TRUE)
ERYOL_dt           = fread("Input/vimc_LE.csv",header = T)  
countries          = fread("Input/vimc_countries.csv",header = T)
age_dist           = fread("Input/vimc_agedist.csv",header = T)
timeliness_pinec   = fread("Input/vimc_timeliness.csv",header = T)

#not country specific 
efficacy_all      = fread("Input/efficacy_all.csv", header = T)
disease_efficacy  = fread("Input/pathogen_disease_efficacy.csv", header = T)  
DALYs_dt          = fread("Input/DALYWEIGHTS.csv",header = T)



# Model set-up
pathogen_list = c('Rota','Sp','Hib')
coverage_type = c("default", "best")

dist_type      = "A1"
psa            = F
waning         = T #currently only Rota Vaccine Efficacy wanes in the model even if you select T for "Sp" and "Hib"
p              = pathogen_list[1]
coverage_given = coverage_type[2]

if("cluster_using_openmpi" %in% commandArgs() & psa == T){
  using_openmpi <- TRUE
  require("doMPI")
}else{
  using_openmpi <- FALSE
}


cov_type = coverage_given
if(coverage_given == coverage_type[1]){
  coverage_dat = fread("Input/corrected_cov_default.csv",header = T)
}else if(coverage_given == coverage_type[2]){
  coverage_dat = fread("Input/corrected_cov_best.csv",header = T)   
}

pathogen_name = p
print(p)

seed.no = which(p == pathogen_list)  #needed for the PSA runs
each_country <- sort(unique(countries$Country))
#each_country <- c("Pakistan", "India", "Nigeria", "Ethiopia", "China")  ## PINE COUNTRIES ONLY
#each_country <- each_country[1]


#cluster set-up
if(using_openmpi == FALSE){
  cores=detectCores()
  use_cluster = 4
  if(use_cluster > 1){
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
  }
  
}else{
  cores=detectCores()
  use_cluster = cores[1]-1
  if(using_openmpi | use_cluster > 1){
    if(!using_openmpi){
      if( use_cluster != detectCores() ){
        warning(paste0(detectCores(), " cores detected but ", use_cluster, " specified."))
      }
      cl <- makeCluster(use_cluster)
      registerDoParallel(cl)
    } else {
      cl <- startMPIcluster()
      registerDoMPI(cl)
      print(paste0("Clustersize: ",clusterSize(cl)))
    }
  }
}


if(psa == T){
  run_until = 200  #Set number of PSA runs here:
  scenario = "PSA"
}else if(psa == F){
  run_until = 1 
  scenario = "Deterministic"
}

writelog("univac_log",paste0("Main; ",scenario, " ", p, " ", coverage_given, " run initiated.")) #Change this to suit both deter or psa runs

foreach(i = each_country, .packages= c('doParallel', 'parallel', 'data.table', 'pbapply','extraDistr')) %dopar% {
  
  country_run = 0
  rand.num <- list()
  set.seed(seed.no)  #used mostly to keep track of PSA runs for each pathogen
  Random.no.list = runif(1000000, 0, 1)
  random.no = 0
  
  stochastic_temp = data.frame()
  rates_dat = data.frame()
  ageDist_df = data.frame()
  timelinessDist_df = data.frame()
  efficacy_PSA = data.frame()
  efficacy_PSA2 = data.frame()
  duration_of_illness = data.frame()
  healthy_time_lost = data.frame()
  duration_of_illness_VAX = data.frame()
  healthy_time_lost_VAX = data.frame()
  .GlobalEnv$i <- i
  
  writelog("univac_log",paste0(scenario, " : Pathogen = ", p,"; ", "Coverage = ",coverage_given))
  
  foreach(run = (1:run_until), .combine = 'c', .packages= c('doParallel', 'parallel', 'data.table', 'pbapply', 
                                                            'extraDistr','Rcpp', .noexport = "mypackage")) %dopar% {
                                                              
                                                              #run = 1   
                                                              .GlobalEnv$run <- run
                                                              country_code = countries[Country == i, country_code]
                                                              assign("country_code", country_code, envir = .GlobalEnv)
                                                              random.no = random.no + 1
                                                              probability = Random.no.list[random.no]
                                                              
                                                              country_run = country_run + 1
                                                              #cat(country_run,".", i, "; Coverage = ",coverage_given, "\n")
                                                              
                                                              if(psa == T){
                                                                mid.pop.country = mid.pop[Country == i]
                                                                low.pop.country = low.pop[Country == i]
                                                                high.pop.country = high.pop[Country == i]
                                                                
                                                                Mid_0 =   mid.pop.country$X0
                                                                Low_0 =   low.pop.country$X0
                                                                High_0 = high.pop.country$X0
                                                                
                                                                Mid_1 =   mid.pop.country$X1
                                                                Low_1 =   low.pop.country$X1
                                                                High_1 = high.pop.country$X1
                                                                
                                                                Mid_2 =   mid.pop.country$X2
                                                                Low_2 =   low.pop.country$X2
                                                                High_2 = high.pop.country$X2
                                                                
                                                                Mid_3 =   mid.pop.country$X3
                                                                Low_3 =   low.pop.country$X3
                                                                High_3 = high.pop.country$X3
                                                                
                                                                Mid_4 =   mid.pop.country$X4
                                                                Low_4 =   low.pop.country$X4
                                                                High_4 = high.pop.country$X4
                                                                
                                                                Mid_Cohort =   mid.pop.country$cohort_pop
                                                                Low_Cohort =   low.pop.country$cohort_pop
                                                                High_Cohort = high.pop.country$cohort_pop
                                                                
                                                                X0 = mapply(Mid_0, Low_0, High_0, probability, FUN = Beta.Pert.dist)
                                                                X1 = mapply(Mid_1, Low_1, High_1, probability, FUN = Beta.Pert.dist)
                                                                X2 = mapply(Mid_2, Low_2, High_2, probability, FUN = Beta.Pert.dist)
                                                                X3 = mapply(Mid_3, Low_3, High_3, probability, FUN = Beta.Pert.dist)
                                                                X4 = mapply(Mid_4, Low_4, High_4, probability, FUN = Beta.Pert.dist)
                                                                cohort_pop = mapply(Mid_Cohort, Low_Cohort, High_Cohort, probability, FUN = Beta.Pert.dist)
                                                                country_pop = as.data.frame(cbind(run.id=run, probability, mid.pop.country[,2:4], X0, X1, X2, X3, X4, cohort_pop))
                                                              }else{
                                                                country_pop = cbind(run.id=run, mid.pop[Country == i])
                                                              }

                                                              country_pop = as.data.table(country_pop)
                                                              POP.2000.PSA = country_pop[Year == 2000]$cohort_pop
                                                              
                                                              #1.
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "country_pop", run.num = run, random.no = random.no, probability = probability))
                                                              
                                                              country = i
                                                              .GlobalEnv$country <- country
                                                              mortality_strata = countries[Country == i, mortality_group]
                                                              #country_pop1 = pop_dt2[Country == i]
                                                              country_rates = rates_all[country == i & pathogen == pathogen_name]
                                                              country_u5mr = u5mr_dat[Country == i]
                                                              country_coverage = coverage_dat[country == i]
                                                              country_timeliness = timeliness_pinec[country == i]
                                                              ERYOL = ERYOL_dt[country == i]
                                                              age_dist_country = age_dist[pathogen == pathogen_name & distribution == dist_type & country == i]
                                                              efficacy_dt = (efficacy_all[mortality_group == mortality_strata])
                                                              
                                                              random.no = random.no + 1
                                                              probability = Random.no.list[random.no]
                                                              rates_df = rates(country_pop, country_rates, country_u5mr, psa, run.id = run, range.type = "Mid", probability)
                                                              rates_dat <- rbind(rates_dat, rates_df)
                                                              
                                                              #2.
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "rates", run.num = run, random.no = random.no, probability = probability))
                                                              
                                                              # NO VACCINATION
                                                              # Age Distribution
                                                              
                                                              # Shape1 and Scale only change for Rota; else only they don't change regardless of Mid, Low and High
                                                              random.no = random.no + 1
                                                              probability = Random.no.list[random.no]
                                                              
                                                              shape1 = round(as.numeric(age_dist_country[param == 'Shape 1' & type == 'Mid', value]), 2)
                                                              shape2 = round(as.numeric(age_dist_country[param == 'Shape 2' & type == 'Mid', value]), 2)
                                                              
                                                              Scale_Mid = round(as.numeric(age_dist_country[param == 'Scale' & type == 'Mid', value]), 2)
                                                              Scale_Low = round(as.numeric(age_dist_country[param == 'Scale' & type == 'Low', value]), 2)
                                                              Scale_High = round(as.numeric(age_dist_country[param == 'Scale' & type == 'High', value]), 2)
                                                              
                                                              if(psa == T){
                                                                Scale = Beta.Pert.dist(Scale_Mid, Scale_Low, Scale_High, probability)
                                                              }else if(psa == F){
                                                                Scale = round(as.numeric(age_dist_country[param == 'Scale' & type == 'Mid', value]), 2)
                                                              }
                                                              
                                                              
                                                              ageDist_df <- rbind(ageDist_df, cbind(country, run, random.no, probability, shape1, shape2, Scale))
                                                              dist = weekly_dist(260, Scale=Scale, shape1=shape1, shape2=shape2)
                                                              
                                                              #3.
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "Age distribution", run.num = run, random.no = random.no, probability = probability))
                                                              
                                                              scenarios = c('Cases', 'Visits', 'Hosps', 'Deaths')
                                                              #Cases before vacciantion
                                                              NoVax_cases_weekly = weekly_scenario_func(rates_df, dist, scenario_type = scenarios[1])
                                                              NoVax_cases_yearly = yearly_scenario_func(rates_df, dist_df = dist, scenario_type = scenarios[1])
                                                              
                                                              #Deaths before vaccination
                                                              NoVax_deaths_weekly = weekly_scenario_func(rates_df, dist, scenario_type = scenarios[4])
                                                              NoVax_deaths_yearly = yearly_scenario_func(rates_df, dist_df = dist, scenario_type = scenarios[4])  #convert weekly deaths to yearly
                                                              
                                                              #Dalys before vaccination - these are the daly's for cohort year 
                                                              random.no = random.no + 1
                                                              probability = Random.no.list[random.no]
                                                              Dalys_PREVAX = as.data.table(t(apply(NoVax_cases_yearly, 1, function(x) Dalys_preVax(x, NoVax_deaths_yearly, ERYOL, DALYs_dt, run.id, random.no, probability, Random.no.list, psa))))
                                                              #4.
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "Dalys PreVax - Healthy time lost", run.num = run, random.no = random.no, probability = probability))
                                                              
                                                              random.no = random.no + 1
                                                              probability = Random.no.list[random.no]
                                                              
                                                              #5.
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "Dalys PreVax - Duration of illness", run.num = run, random.no = random.no, probability = probability))
                                                              
                                                              Dalys_PreV = t(Dalys_PREVAX)
                                                              
                                                              dalys_PSA1_output <- list()
                                                              list_index = seq(1, length(Dalys_PreV), 39)
                                                              for(j in 1:length(list_index)){
                                                                index = list_index[j]
                                                                dalys_PSA1_output[j] <- Dalys_PreV[[index]]["healthy_time_lost_PSA"]
                                                              }
                                                              dalys_PSA1_output =  do.call("rbind", dalys_PSA1_output)
                                                              healthy_time_lost = rbind(healthy_time_lost, dalys_PSA1_output)
                                                              
                                                              dalys_PSA2_output <- list()
                                                              list_index = seq(1, length(Dalys_PreV), 39)
                                                              for(j in 1:length(list_index)){
                                                                index = list_index[j]
                                                                dalys_PSA2_output[j] <- Dalys_PreV[[index]]["duration_of_illness_PSA"]
                                                              }
                                                              dalys_PSA2_output =  do.call("rbind", dalys_PSA2_output)
                                                              duration_of_illness = rbind(duration_of_illness, dalys_PSA2_output)
                                                              
                                                              Dalys_PREVAX <- list()
                                                              for(d in 1:length(Dalys_PreV)){
                                                                Dalys_PREVAX[d] <- Dalys_PreV[[d]]["dalys_preVax"]
                                                              }
                                                              Dalys_PREVAX =  as.data.table(do.call("rbind", Dalys_PREVAX))
                                                              
                                                              Dalys_PREVAX = Dalys_PREVAX[, (names(Dalys_PREVAX[,3:8])) := lapply(.SD, as.numeric), .SDcols = names(Dalys_PREVAX[,3:8])]  #coverting all the numbers to numeric values
                                                              Dalys_PreVAC = Dalys_PREVAX[, lapply(.SD, sum), by = c(names(Dalys_PREVAX[,3])), .SDcols = (names(Dalys_PREVAX[,4:8]))]
                                                              Dalys_PreVAC = cbind(country = country_code, pathogen= DALYs_dt[disease == as.character(NoVax_cases_yearly[1,1]), pathogen], Dalys_PreVAC)
                                                              
                                                              #PSA - timeliness  
                                                              random.no = random.no + 1
                                                              probability = Random.no.list[random.no]
                                                              
                                                              
                                                              if(psa == T){
                                                                shape_d1_mid = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Shape' & type == 'Mid', value])
                                                                shape_d1_low = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Shape' & type == 'Low', value])
                                                                shape_d1_high = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Shape' & type == 'High', value])
                                                                shape_d1 = Beta.Pert.dist(shape_d1_mid, shape_d1_low, shape_d1_high, probability)
                                                                
                                                                Scale_d1_Mid = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Scale' & type == 'Mid', value])
                                                                Scale_d1_Low = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Scale' & type == 'Low', value])
                                                                Scale_d1_High = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Scale' & type == 'High', value])
                                                                Scale_d1 = Beta.Pert.dist(Scale_d1_Mid, Scale_d1_Low, Scale_d1_High, probability)
                                                                
                                                                shape_d2_mid = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Shape' & type == 'Mid', value])
                                                                shape_d2_low = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Shape' & type == 'Low', value])
                                                                shape_d2_high = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Shape' & type == 'High', value])
                                                                shape_d2 = Beta.Pert.dist(shape_d2_mid, shape_d2_low, shape_d2_high, probability)
                                                                
                                                                Scale_d2_Mid = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Scale' & type == 'Mid', value])
                                                                Scale_d2_Low = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Scale' & type == 'Low', value])
                                                                Scale_d2_High = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Scale' & type == 'High', value])
                                                                Scale_d2 = Beta.Pert.dist(Scale_d2_Mid, Scale_d2_Low, Scale_d2_High, probability)
                                                                
                                                                shape_d3_mid = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Shape' & type == 'Mid', value])
                                                                shape_d3_low = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Shape' & type == 'Low', value])
                                                                shape_d3_high = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Shape' & type == 'High', value])
                                                                shape_d3 = Beta.Pert.dist(shape_d3_mid, shape_d3_low, shape_d3_high, probability)
                                                                
                                                                Scale_d3_Mid = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Scale' & type == 'Mid', value])
                                                                Scale_d3_Low = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Scale' & type == 'Low', value])
                                                                Scale_d3_High = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Scale' & type == 'High', value])
                                                                Scale_d3 = Beta.Pert.dist(Scale_d3_Mid, Scale_d3_Low, Scale_d3_High, probability)
                                                                
                                                              }else if(psa == F){
                                                                shape_d1 = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Shape' & type == 'Mid', value])
                                                                Scale_d1 = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Scale' & type == 'Mid', value])
                                                                
                                                                shape_d2 = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Shape' & type == 'Mid', value])
                                                                Scale_d2 = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Scale' & type == 'Mid', value])
                                                                
                                                                shape_d3 = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Shape' & type == 'Mid', value])
                                                                Scale_d3 = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Scale' & type == 'Mid', value])
                                                              }
                                                              
                                                              Shift_d1 = as.numeric(country_timeliness[dose == 'DTP1' & param == 'Shift' & type == 'Mid', value])
                                                              Shift_d2 = as.numeric(country_timeliness[dose == 'DTP2' & param == 'Shift' & type == 'Mid', value])
                                                              Shift_d3 = as.numeric(country_timeliness[dose == 'DTP3' & param == 'Shift' & type == 'Mid', value])

                                                              
                                                              d1_timeliness = start_weekly_dist_shift(260, Scale = Scale_d1, shape1 = shape_d1, Shift = Shift_d1)
                                                              d2_timeliness = start_weekly_dist_shift(260, Scale = Scale_d2, shape1 = shape_d2, Shift = Shift_d2)
                                                              d3_timeliness = start_weekly_dist_shift(260, Scale = Scale_d3, shape1 = shape_d3, Shift = Shift_d3)
                                                              
                                                              timelinessDist_df <- rbind(timelinessDist_df, cbind(country, run, random.no, probability, shape_d1, Shift_d1, Scale_d1, shape_d2, Shift_d2, Scale_d2, shape_d3, Shift_d3, Scale_d3))
                                                              #6.
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "timeliness", run.num = run, random.no = random.no, probability = probability))
                                                              
                                                              
                                                              #7.
                                                              random.no = random.no + 1
                                                              probability = Random.no.list[random.no]
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "Efficacy d1/d2/d3", run.num = run, random.no = random.no, probability = probability))
                                                              
                                                              Cases_PreVax_dat  = NoVax_cases_weekly[ ! NoVax_cases_weekly$Y %in% c(1996:1999), ]
                                                              impact_table <- data.table(matrix(0, nrow = 1, ncol = 262))  #create impact table to store the impact calculations by E1/E2 type and year
                                                              
                                                              writelog("univac_log",paste0(scenario, " : ", country_code, "; Run ",run,"/",run_until,"; Start ", p, " ",coverage_given, " model"))
                                                              lst_output = as.data.table(t(pbapply(as.matrix(Cases_PreVax_dat), 1, function(x) Vax_impact(pathogen_name, run.id = run, random.no, probability, Random.no.list, country_code,
                                                                                                                                                          country, x, NoVax_deaths_weekly,country_pop,vaccine_start=2000,
                                                                                                                                                          country_coverage,d1_timeliness, d2_timeliness, d3_timeliness,
                                                                                                                                                          disease_efficacy, efficacy_dt, ERYOL, DALYs_dt,impact_table,psa,waning))))
                                                              #8/9/10.
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "Alpha/Mean", run.num = run, random.no = (random.no+1), probability = Random.no.list[random.no+1]))
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "Dalys PostVax - Healthy time lost", run.num = run, random.no = (random.no+2), probability = Random.no.list[random.no+2]))
                                                              rand.num <- rbind(rand.num, c(country.code = country_code, variable = "Dalys PostVax - Duration of illness", run.num = run, random.no = (random.no+3), probability = Random.no.list[random.no+3]))
                                                              
                                                              
                                                              random.no = random.no + 3
                                                              probability = Random.no.list[random.no]
                                                              
                                                              lst_out = t(lst_output)
                                                              
                                                              if(pathogen_name == "Hib"){
                                                                list_n = 106
                                                              }else{
                                                                list_n = 36
                                                              }
                                                              
                                                              efficacy_PSA_output <- list()
                                                              list_index = c(1,list_n)
                                                              for(j in 1:length(list_index)){
                                                                index = list_index[j]
                                                                efficacy_PSA_output[j] <- lst_out[[index]]["efficacy_PSA_output"]
                                                              }
                                                              efficacy_PSA_output =  do.call("rbind", efficacy_PSA_output)
                                                              efficacy_PSA = rbind(efficacy_PSA, efficacy_PSA_output)
                                                              
                                                              
                                                              efficacyPSA2 <- list()
                                                              list_index = c(1,list_n)
                                                              for(j in 1:length(list_index)){
                                                                index = list_index[j]
                                                                efficacyPSA2[j] <- lst_out[[index]]["efficacyPSA2"]
                                                              }
                                                              efficacyPSA2 =  do.call("rbind", efficacyPSA2)
                                                              efficacy_PSA2 = rbind(efficacy_PSA2, efficacyPSA2)
                                                              
                                                              
                                                              Dalys_PSA1_VAX <- list()
                                                              list_index = seq(1, length(lst_out), 35)
                                                              for(j in 1:length(list_index)){
                                                                index = list_index[j]
                                                                Dalys_PSA1_VAX[j] <- lst_out[[index]]["healthy_time_lost_PSA"]
                                                              }
                                                              Dalys_PSA1_VAX =  do.call("rbind", Dalys_PSA1_VAX)
                                                              healthy_time_lost_VAX = rbind(healthy_time_lost_VAX, Dalys_PSA1_VAX)
                                                              
                                                              
                                                              Dalys_PSA2_VAX <- list()
                                                              list_index = seq(1, length(lst_out), 35)
                                                              for(j in 1:length(list_index)){
                                                                index = list_index[j]
                                                                Dalys_PSA2_VAX[j] <- lst_out[[index]]["duration_of_illness_PSA"]
                                                              }
                                                              Dalys_PSA2_VAX =  do.call("rbind", Dalys_PSA2_VAX)
                                                              duration_of_illness_VAX = rbind(duration_of_illness_VAX, Dalys_PSA2_VAX)
                                                              
                                                              
                                                              suppressPackageStartupMessages(library(tidyverse))
                                                              #https://community.rstudio.com/t/spread-with-multiple-value-columns/5378
                                                              
                                                              rates_2000 = rates_df[Y == 2000]
                                                              rates_2000 = rates_2000[,-c(3,4,6)]
                                                              
                                                              rates_spread <- rates_2000 %>%
                                                                gather(variable, value, -(country:Disease)) %>%
                                                                unite(temp, Disease, variable) %>%
                                                                spread(temp, value)
                                                              efficacy_spread = as.data.frame(efficacy_PSA_output[,-c(1:5)])
                                                              efficacy_spread <- efficacy_spread %>%
                                                                gather(variable, value, -(efficacy_type)) %>%
                                                                unite(temp, efficacy_type, variable) %>%
                                                                spread(temp, value)
                                                              
                                                              efficacy_spread2 = as.data.frame(efficacyPSA2[,-c(1:5)])
                                                              efficacy_spread2 <- efficacy_spread2 %>%
                                                                gather(variable, value, -(efficacy_type)) %>%
                                                                unite(temp, efficacy_type, variable) %>%
                                                                spread(temp, value)
                                                              
                                                              Dalys_spread = as.data.frame(dalys_PSA1_output)[,-c(1:5)]
                                                              Dalys_spread <- Dalys_spread %>%
                                                                gather(variable, value, -(disease_type)) %>%
                                                                unite(temp, disease_type, variable) %>%
                                                                spread(temp, value)
                                                              Dalys_spread2 = as.data.frame(dalys_PSA2_output)[,-c(1:5)]
                                                              Dalys_spread2 <- Dalys_spread2 %>%
                                                                gather(variable, value, -(disease_type)) %>%
                                                                unite(temp, disease_type, variable) %>%
                                                                spread(temp, value)
                                                              Dalys_vax_spread = as.data.frame(Dalys_PSA1_VAX)[,-c(1:5)]
                                                              Dalys_vax_spread <- Dalys_vax_spread %>%
                                                                gather(variable, value, -(disease_type)) %>%
                                                                unite(temp, disease_type, variable) %>%
                                                                spread(temp, value)
                                                              Dalys_vax_spread2 = as.data.frame(Dalys_PSA2_VAX)[,-c(1:5)]
                                                              Dalys_vax_spread2 <- Dalys_vax_spread2 %>%
                                                                gather(variable, value, -(disease_type)) %>%
                                                                unite(temp, disease_type, variable) %>%
                                                                spread(temp, value)
                                                              
                                                              #find a way to output efficacy in this format so that it's not storing excess data & keep consistent prevax and postvax outcome names so that it's not storing multiple data frames.
                                                              
                                                              if(psa == T){
                                                                stochastic_temp =  rbind(stochastic_temp, cbind(rates_spread, POP.2000.PSA, ageDist_shape1=shape1, ageDist_shape2=shape2, ageDist_scale=Scale,
                                                                                                                timeliness_d1_shape=shape_d1, timeliness_d1_shift=Shift_d1, timeliness_d1_scale=Scale_d1,
                                                                                                                timeliness_d2_shape=shape_d2, timeliness_d2_shift=Shift_d2, timeliness_d2_scale=Scale_d2,
                                                                                                                timeliness_d3_shape=shape_d3, timeliness_d3_shift=Shift_d3, timeliness_d3_scale=Scale_d3,
                                                                                                                preVax=Dalys_spread, preVax=Dalys_spread2, Dalys_vax_spread, Dalys_vax_spread2, 
                                                                                                                efficacy_spread, efficacy_spread2))
                                                                
                                                                assign("stochastic_temp", stochastic_temp, envir = .GlobalEnv)
                                                              }
                                                              
                                                              # --------------------------------  need to edit from here to save each output from PSA for each country  --------------------------------------
                                                              
                                                              all_scenarios <- list()
                                                              for(n in 1:length(lst_out)){
                                                                all_scenarios[n] <- lst_out[[n]]["combined_output"]
                                                              }
                                                              all_scenarios = do.call("rbind", all_scenarios)
                                                              
                                                              # This included the dalys for cohort year
                                                              all_scenarios = as.data.table(all_scenarios)
                                                              all_scenarios = all_scenarios[, (names(all_scenarios[,5:24])) := lapply(.SD, as.numeric), .SDcols = names(all_scenarios[,5:24])]  #coverting all the numbers to numeric values
                                                              post_vaccination = all_scenarios[, lapply(.SD, sum), by = c(names(all_scenarios[,2:9])), .SDcols = (names(all_scenarios[,10:24]))]
                                                              post_vaccination = cbind(pathogen = pathogen_name, post_vaccination)
                                                              
                                                              
                                                              # Cases_PreVaxx = NoVax_cases_yearly[  NoVax_cases_yearly$Y %in% c(1996:1999), ]
                                                              # Deaths_PreVaxx = NoVax_deaths_yearly[  NoVax_deaths_yearly$Y %in% c(1996:1999), ]
                                                              
                                                              PreVaxx1 = cbind(NoVax_cases_yearly, NoVax_deaths_yearly)
                                                              PreVaxx1 = PreVaxx1[, -c(3,9,10,11)]
                                                              PreVaxx1 = as.data.table(PreVaxx1)
                                                              PreVaxx_cases_deaths1 = PreVaxx1[, lapply(.SD, sum), by = c(names(PreVaxx1[,2])), .SDcols = (names(PreVaxx1[,3:12]))]
                                                              #Dalys_PreVaxx =  Dalys_PreVAC[  Dalys_PreVAC$year %in% c(1996:1999), ]
                                                              pre_vaccination = cbind(PreVaxx_cases_deaths1, Dalys_PreVAC)
                                                              pre_vaccination =  pre_vaccination[,c(-1)]
                                                              pre_vaccination = cbind(pre_vaccination, country)
                                                              colnames(pre_vaccination)[1:10] = c("cases_0", "cases_1","cases_2","cases_3","cases_4","deaths_0","deaths_1","deaths_2","deaths_3","deaths_4")
                                                              colnames(pre_vaccination)[19] = c("country_name")
                                                              cohort_preVax = country_pop[Year %in% c(1996:2034)][,6:10]
                                                              colnames(cohort_preVax) = c(colnames(post_vaccination)[5:9])
                                                              pre_vaccination = cbind(cohort_preVax, pre_vaccination)
                                                              pre_vaccination = setcolorder(pre_vaccination, (colnames(post_vaccination)))
                                                              
                                                              #pre_vaccination = all_scenarios_PreVax[ ! all_scenarios_PreVax$year %in% c(2031:2032), ]
                                                              post_vaccination = rbind(pre_vaccination[year %in% c(1996:1999)], post_vaccination)
                                                              
                                                              pre_vaccination = diagonal_pop(pre_vaccination)
                                                              post_vaccination = diagonal_pop(post_vaccination)
                                                              
                                                              
                                                              # #transposing table from wide to long format
                                                              #Pre Vaccination
                                                              prevax_output = reshape(data = pre_vaccination,
                                                                                      idvar = c("pathogen","year","country"),
                                                                                      timevar = 'age',
                                                                                      times = c(0:4),
                                                                                      varying = list(cohort_size=c(5:9),cases=c(10:14),deaths=c(15:19),dalys=c(20:24)),
                                                                                      direction="long",
                                                                                      v.names = c("cohort_size","cases","deaths","dalys"),
                                                                                      sep="_")
                                                              
                                                              colnames(prevax_output)[1] <- "disease"
                                                              prevax_output <- setcolorder(prevax_output, c("disease", "year", "age", "country", "country_name", "cohort_size", "cases", "dalys", "deaths"))
                                                              
                                                              prevax_output = prevax_output[!(prevax_output$age == "0" & prevax_output$year > 2030 ),]
                                                              prevax_output = prevax_output[!(prevax_output$age == "1" & prevax_output$year > 2031 ),]
                                                              prevax_output = prevax_output[!(prevax_output$age == "2" & prevax_output$year > 2032 ),]
                                                              prevax_output = prevax_output[!(prevax_output$age == "3" & prevax_output$year > 2033 ),]
                                                              
                                                              prevax_output = cbind(run_id=run, prevax_output)
                                                          
                                                              
                                                              #Post Vaccination
                                                              desired_montague_template = reshape(data = post_vaccination,
                                                                                                  idvar = c("pathogen","year","country"),
                                                                                                  timevar = 'age',
                                                                                                  times = c(0:4),
                                                                                                  varying = list(cohort_size=c(5:9),cases=c(10:14),deaths=c(15:19),dalys=c(20:24)),
                                                                                                  direction="long",
                                                                                                  v.names = c("cohort_size","cases","deaths","dalys"),
                                                                                                  sep="_")
                                                              
                                                              colnames(desired_montague_template)[1] <- "disease"
                                                              desired_montague_template <- setcolorder(desired_montague_template, c("disease", "year", "age", "country", "country_name", "cohort_size", "cases", "dalys", "deaths"))
                                                              
                                                              desired_montague_template = desired_montague_template[!(desired_montague_template$age == "0" & desired_montague_template$year > 2030 ),]
                                                              desired_montague_template = desired_montague_template[!(desired_montague_template$age == "1" & desired_montague_template$year > 2031 ),]
                                                              desired_montague_template = desired_montague_template[!(desired_montague_template$age == "2" & desired_montague_template$year > 2032 ),]
                                                              desired_montague_template = desired_montague_template[!(desired_montague_template$age == "3" & desired_montague_template$year > 2033 ),]
                                                              
                                                              desired_montague_template = cbind(run_id=run, desired_montague_template)
                                                              
                                                              prevax_output[, 8:10] <- round(prevax_output[, 8:10], digits = 3)
                                                              desired_montague_template[, 8:10] <- round(desired_montague_template[, 8:10], digits = 3)
                                                              
                                                              assign("prevax_output", prevax_output, envir = .GlobalEnv)
                                                              assign("desired_montague_template", desired_montague_template, envir = .GlobalEnv)
                                                              assign("rand.num", rand.num, envir = .GlobalEnv)
                                                              
                                                              if(psa == T){
                                                                if(!dir.exists(paste0("Output/VIMC_", pathogen_name, "_Results/","PSA"))){
                                                                  dir.create(paste0("Output/VIMC_", pathogen_name, "_Results/","PSA"))
                                                                }
                                                                
                                                                if(!dir.exists(paste0("Output/VIMC_", pathogen_name, "_Results/PSA/",pathogen_name, "_", Sys.Date()))){
                                                                  dir.create(paste0("Output/VIMC_", pathogen_name, "_Results/PSA/",pathogen_name, "_", Sys.Date()))
                                                                }
                                                                path_to_save = (paste0("Output/VIMC_", pathogen_name, "_Results/PSA/",pathogen_name, "_", Sys.Date(),"/"))
                                                                assign("path_to_save", path_to_save, envir = .GlobalEnv)
                                                                
                                                                
                                                                if(coverage_given == coverage_type[1]){
                                                                  save_novax = dir.create(paste0(path_to_save,"no_vax"))
                                                                  if(!dir.exists(paste0(path_to_save,"no_vax/run", run))){
                                                                    dir.create(paste0(path_to_save,"no_vax/run", run))
                                                                  }
                                                                  path_to_save_no_vax = (paste0(path_to_save,"no_vax/run", run, "/"))
                                                                  
                                                                  
                                                                  save_scenario = dir.create(paste0(path_to_save,"default_routine"))
                                                                  if(!dir.exists(paste0(path_to_save,"default_routine/run", run))){
                                                                    dir.create(paste0(path_to_save,"default_routine/run", run))
                                                                  }
                                                                  path_to_save_scenario = (paste0(path_to_save,"default_routine/run", run, "/"))
                                                                  fwrite(prevax_output, file = paste0(path_to_save_no_vax, paste0(country_code, "_", pathogen_name,"_", "NoVax_run", run,".csv")), sep = ",", col.names = TRUE, row.names = FALSE)
                                                                  
                                                                }else{
                                                                  
                                                                  dir.create(paste0(path_to_save,"best_routine"))
                                                                  if(!dir.create(paste0(path_to_save,"best_routine/run", run))){
                                                                    dir.create(paste0(path_to_save,"best_routine/run", run))
                                                                  }
                                                                  path_to_save_scenario = (paste0(path_to_save,"best_routine/run", run, "/"))
                                                                }
                                                                
                                                                fwrite(desired_montague_template, file = paste0(path_to_save_scenario, paste0(country_code, "_", pathogen_name,"_", cov_type, "_run",run,".csv")), sep = ",", col.names = TRUE, row.names = FALSE)
                                                                stochastic_temp = as.data.table(stochastic_temp)
                                                                fwrite(stochastic_temp, file = paste0(path_to_save, paste0(country_code, "_PSA_", pathogen_name,"_InputParams",".csv")), sep = ",", col.names = TRUE, row.names = FALSE)
                                                                fwrite(rand.num, file = paste0(path_to_save, paste0(country_code, "_PSA_", pathogen_name,"_RandomNumbers",".csv")), sep = ",", col.names = TRUE, row.names = FALSE) 
                                                                
                                                              }else if(psa == F){
                                                                
                                                                if(!dir.exists(paste0("Output/VIMC_", pathogen_name, "_Results/","deter"))){
                                                                  dir.create(paste0("Output/VIMC_", pathogen_name, "_Results/","deter"))
                                                                }
                                                                
                                                                if(!dir.exists(paste0("Output/VIMC_", pathogen_name, "_Results/deter/",pathogen_name, "_deter_", Sys.Date()))){
                                                                  dir.create(paste0("Output/VIMC_", pathogen_name, "_Results/deter/",pathogen_name, "_deter_", Sys.Date()))
                                                                }
                                                                path_to_save = (paste0("Output/VIMC_", pathogen_name, "_Results/deter/",pathogen_name, "_deter_", Sys.Date(),"/"))
                                                                
                                                                
                                                                if(coverage_given == coverage_type[1]){
                                                                  
                                                                  save_novax = dir.create(paste0(path_to_save,"no_vax"))
                                                                  path_to_save_no_vax = (paste0(path_to_save,"no_vax/"))

                                                                  save_scenario = dir.create(paste0(path_to_save,"default_routine"))
                                                                  path_to_save_scenario = (paste0(path_to_save,"default_routine/"))
                                                                  
                                                                  if(pathogen_name == "Rota"){
                                                                    fwrite(prevax_output, file = paste0(path_to_save_no_vax, paste0(country_code, "_",pathogen_name,"_", "NoVax_deter_waning", waning, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE)
                                                                  }else{
                                                                    fwrite(prevax_output, file = paste0(path_to_save_no_vax, paste0(country_code, "_",pathogen_name,"_", "NoVax_deter.csv")), sep = ",", col.names = TRUE, row.names = FALSE)
                                                                  }
                                                                  
                                                                }else{
                                                                  dir.create(paste0(path_to_save,"best_routine"))
                                                                  path_to_save_scenario = (paste0(path_to_save,"best_routine/"))
                                                                }
                                                                
                                                                if(pathogen_name == "Rota"){
                                                                  fwrite(desired_montague_template, file = paste0(path_to_save_scenario, paste0(country_code, "_",pathogen_name,"_", cov_type,"_deter_waning", waning, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE)
                                                                }else{
                                                                  fwrite(desired_montague_template, file = paste0(path_to_save_scenario, paste0(country_code, "_",pathogen_name,"_", cov_type,"_deter.csv")), sep = ",", col.names = TRUE, row.names = FALSE)
                                                                }

                                                                }
                                                                
                                                                writelog("univac_log",paste0(scenario, " : ", country_code, "; Run ",run,"/",run_until,"; Finished ", p, " ",coverage_given, " model"))
                                                              }

                                                            }


#clean environment
writelog("univac_log",paste0(scenario, " ", p, " ", coverage_given, " coverage finished."))
if(using_openmpi | use_cluster > 1){
  stopCluster(cl)
}
if(using_openmpi){
  mpi.quit()
}


