#rm(list=ls())
library(data.table)
library(dplyr)
library(tictoc)
library(Rcpp)

#1. small function to write message a to logfile logname
writelog <- function(logname, x){
  write(
    paste0(
      format(Sys.time(), "%Y/%m/%d %H:%M:%S"),
      " ",
      x
    ),
    file=logname,
    append=TRUE
  )
}


#2. Function returns the sum of populations <5 for each cohort year (gets diagonal sum)
cohort_pop <- function(df) {
  library(data.table)
  print(df)
  pop_all = data.table()
  
  for(i in unique(df$Country)){
    #subset pop data for each country
    country_subset = df[Country == i]  
    
    #convert pop df to matrix excluding year 
    pop_matrix <- (as.matrix(country_subset[, -(1:5)]))
    
    #gets the diagonal indicator of the matrix
    pop_indicator <- row(pop_matrix) - col(pop_matrix)
    
    #gets a list of lists containing populations changes for each cohort year
    cohort_year <- split(pop_matrix, pop_indicator)
    
    #gets the total number of individuals in each cohort year
    sum_cohort <- data.frame(unlist(lapply(cohort_year, function(x) sum(x))))
    
    #discards the sum of cohorts before 1996 (anyhting sum before index 0)
    sum_cohort1 <- subset(sum_cohort, as.numeric(rownames(sum_cohort))>=0)
    
    #changing the column name
    colnames(sum_cohort1) <- c("cohort_pop")
    
    #reset row index
    rownames(sum_cohort1)=NULL
    
    pop_df <- as.data.table(cbind(country_subset, sum_cohort1))
    pop_all = rbind(pop_all, pop_df)
  }
  
  return(pop_all)
  
}

#3. Gives the diagonal population matrix
diagonal_pop <- function(df) {
  library(data.table)
  #print(df)
  
  cohort = (as.matrix(df[,c(5:9)]))
  cases = (as.matrix(df[,c(10:14)]))
  deaths = (as.matrix(df[,c(15:19)]))
  dalys = (as.matrix(df[,c(20:24)]))
  
  scenario_lst = list(cohort, cases, deaths, dalys)
  scenario_names = c("cohort", "cases", "deaths", "dalys")
  join_scenarios = list()
  k = 0
  
  for(i in 1:length(scenario_lst)){
    k= k+1
    scenario = scenario_names[i]
    pop_dt = as.data.table(scenario_lst[i])
    pop_matrix = as.matrix(pop_dt)
    
    pop_matrix <- apply(pop_matrix, 2, rev)
    
    #gets the diagonal indicator of the matrix
    pop_indicator <- row(pop_matrix) - col(pop_matrix)
    
    #gets a list of lists containing populations changes for each cohort year
    cohort_year <- split(pop_matrix, pop_indicator)
    
    cohort_scenario = do.call(rbind, cohort_year)
    cohort_scenario = cohort_scenario[c(5:39),]
    
    cohort_scenario_inverse = as.data.table(apply(cohort_scenario, 2, rev))
    colnames(cohort_scenario_inverse) <- colnames(pop_dt)
    
    join_scenarios[[k]] <- cohort_scenario_inverse
  }
  
  big_data = do.call(cbind, join_scenarios)
  data_output = data_output = cbind(df[5:39,1:4], big_data)
  return(data_output)
  
}

#4. Beta Pert distribution
Beta.Pert.dist <- function(Mode, Min.val, Max.val, Rand.probability){
  
  #if(Mode.val == Min.val & Mode.val == Max.val){
  if(isTRUE(all.equal(Mode, Min.val, tolerance = 0.1)) & isTRUE(all.equal(Mode, Max.val, tolerance = 0.1))){
    value = Mode
  }else{
    alpha.val = ((2*(Max.val+4*Mode-5*Min.val))/(3*(Max.val-Min.val)))*(1+4*(((Mode-Min.val)*(Max.val-Mode))/(Max.val-Min.val)^2))
    beta.val  = ((2*(5*Max.val-4*Mode-Min.val))/(3*(Max.val-Min.val)))*(1+4*(((Mode-Min.val)*(Max.val-Mode))/(Max.val-Min.val)^2))
    
    value = round(qnsbeta(Rand.probability, alpha.val, beta.val, min = Min.val, max = Max.val), 2) 
  }
  return(value)
  
}


#5. Specify range.type i.e. "Mid", "Low", "High" for when psa is False
# psa = TRUE will sample from the beta-pert distribution; FALSE will take range.type i.e. Mid range values
# Function returns a data.table of the proportions of cases, visits, hosps and deaths for each disease and cohort year 

rates <- function(country_dt, disease_rates, u5mr_dat, psa, run.id, range.type, probability){
  
  library(data.table)
  disease_rates = as.data.table(disease_rates)
  diseases = unique(disease_rates$disease)  #change column Disease to disease in india rates dt
  years = unique(country_dt$Year)
  
  #Assign empty dataframe and empty lists to store the output from the loop
  rates_dt <- data.frame()
  Cases <- list()
  Visits <- list()
  Hosps <- list()
  Deaths <- list()
  Disease <- list()
  Y <- list()           #empty list to record the year for which each we calculate the proportion of pop under different different disease rates 
  k=0                   #acts as a counter when incrementing years in the loop
  
  #Loop gets the proportion of cases, visits, hosps and deaths using the population for each cohort year
  #there are 7 disease types and 36 years (1996-2031), so there should be 7*36=252 rates for each scenario 
  
  for(j in 1:length(years)){
    
    per_year_df <- country_dt[Year == years[j],]
    pop <- per_year_df$cohort_pop
    u5mr <- u5mr_dat[Year == years[j]]$`%_decrease_from_Ref_Year`
    
    for(i in 1:length(diseases)){
      k = k+1
      
      Y[k] <- years[j]
      
      Disease[i] <- as.character(diseases[i])
      
      case_rates <- disease_rates[disease == diseases[i] & scenario =="Cases",]  #changed from Scenarios to scenario
      visit_rates <- disease_rates[disease == diseases[i] & scenario =="Visits",]
      hosp_rates <- disease_rates[disease == diseases[i] & scenario =="Hosps",]
      death_rates <- disease_rates[disease == diseases[i] & scenario =="Deaths",]
      
      if(psa == FALSE){
        #Disease event rates -- the total number of cases per 100,000 = total <5 pop * rate(n/100,000)
        #Cases[k] <- pop * case_rates$Mid/100000
        #Visits[k] <- pop * visit_rates$Mid/100000
        Cases[k] <- pop * as.numeric(case_rates[, ..range.type])/100000
        Visits[k] <- pop * as.numeric(visit_rates[, ..range.type])/100000
        
        
        if(is.data.frame(hosp_rates) && nrow(hosp_rates)==0){
          Hosps[k] <- 0
        }else{
          #Hosps[k] <- pop * hosp_rates$Mid/100000
          Hosps[k] <- pop * as.numeric(hosp_rates[, ..range.type])/100000
        }
        
        if(is.data.frame(death_rates) && nrow(death_rates)==0){
          Deaths[k] <- 0
        }else{
          #DR <- death_rates$Mid - death_rates$Mid*u5mr   #adjusting the death rate per year
          DR <- death_rates[, ..range.type] - (death_rates[, ..range.type])*u5mr
          Deaths[k] <- pop * DR/100000
        }
        Rate <- range.type
        
      }else if(psa == TRUE){
        
        unscaled_rates_disease = c("Sp pneumonia (non-severe)", "Hib pneumonia (non-severe)")
        
        if(Disease[i] %in% unscaled_rates_disease){
          severe.disease = paste0(stringr::word(Disease[i], 1, 2), " (severe)")
          severe.cases.rates = disease_rates[disease == severe.disease & scenario =="Cases",] 
          severe.visit.rates = disease_rates[disease == severe.disease & scenario =="Visits",]
          
          #adjust rates where mid < low for low and high rates for those particular disease types only
          new.low.case_rates = severe.cases.rates$Low/severe.cases.rates$Mid*case_rates$Mid
          new.high.case_rates = severe.cases.rates$High/severe.cases.rates$Mid*case_rates$Mid
          
          new.low.visit_rates = severe.visit.rates$Low/severe.visit.rates$Mid*visit_rates$Mid
          new.high.visit_rates = severe.visit.rates$High/severe.visit.rates$Mid*visit_rates$Mid
          
          cases_psa_rate = Beta.Pert.dist(Mode = case_rates$Mid, Min.val = new.low.case_rates, Max.val = new.high.case_rates, probability)
          Cases[k] <- pop * cases_psa_rate/100000
          
          visits_psa_rate = Beta.Pert.dist(Mode = visit_rates$Mid, Min.val = new.low.visit_rates, Max.val = new.high.visit_rates, probability)
          Visits[k] <- pop * visits_psa_rate/100000
          
        }else{
          cases_psa_rate = Beta.Pert.dist(Mode = case_rates$Mid, Min.val = case_rates$Low, Max.val = case_rates$High, probability)
          Cases[k] <- pop * cases_psa_rate/100000
          
          visits_psa_rate = Beta.Pert.dist(Mode = visit_rates$Mid, Min.val = visit_rates$Low, Max.val = visit_rates$High, probability)
          Visits[k] <- pop * visits_psa_rate/100000
        }
        
        if(is.data.frame(hosp_rates) && nrow(hosp_rates)==0){
          Hosps[k] <- 0
        } else{
          hosp_psa_rate = Beta.Pert.dist(Mode = hosp_rates$Mid, Min.val = hosp_rates$Low, Max.val = hosp_rates$High, probability)
          Hosps[k] <- pop * hosp_psa_rate/100000
        }
        
        if(is.data.frame(death_rates) && nrow(death_rates)==0){
          Deaths[k] <- 0
        } else{
          DR = Beta.Pert.dist(Mode = death_rates$Mid, Min.val = death_rates$Low, Max.val = death_rates$High, probability) -
            Beta.Pert.dist(Mode = death_rates$Mid, Min.val = death_rates$Low, Max.val = death_rates$High, probability)*u5mr
          
          Deaths[k] <- pop * DR/100000
        }
        Rate <- probability
        
      }
    }
    
  }
  country = country_dt$Country[1]
  rates_dt <- cbind(country, run.id, random.no, Y, Disease, Rate, Cases, Visits, Hosps, Deaths)
  rdt <- as.data.table(rates_dt)
  return(rdt)
}

#6. Burr distribution
burr_func <- function(weekly_age, Scale, shape1, shape2){
  
  burr_dist = 1-((1+(weekly_age/Scale)^shape1)^-shape2)
  
  return(burr_dist)
}


#7. Distribution of cases etc by week of age
weekly_dist <- function(week, Scale, shape1, shape2){
  
  # R package required to run function
  library(data.table)
  
  # 1. CDF end week
  cdf_list <- list()
  Week <- list()
  for(i in 1:week){
    weekly_age = i
    cdf_list[i] = burr_func(weekly_age, Scale, shape1, shape2)
    Week[i] = i
  }
  
  
  # 2. Rescaled PDF
  #Re-scale cdf so that the maximum cumulative frequency adds up to 1
  cdf <- unlist(cdf_list)
  rescaled_list <- list()
  max_cdf = cdf[week]
  
  #rescale cdf by dividing cdf by maximum cdf : cdf_xi/cdf_xmax
  for(i in 1:length(cdf)){
    rescaled_list[i] = cdf[i]/max_cdf
  }
  
  
  # 3. PDF
  # Probability density function by subtracting pdf_x2 = cdf_x2 - cdf_x1
  pdf_list <- list()
  rescaled_cdf <- unlist(rescaled_list)
  pdf_list[1] <- rescaled_cdf[1]    #Keeps the first value of the rescaled cdf
  
  for(i in 2:length(rescaled_cdf)){
    pdf_list[i] <- rescaled_cdf[i] - rescaled_cdf[i-1]
  }
  
  dist_dt <- as.data.table(cbind(Week, CDF=cdf_list, RESCALED_CDF=rescaled_list, PDF=pdf_list))
  return(dist_dt)
}



#8. timeliness in univac is calculated by an additional shift parameter instead of shape
burr_shift_pdf <- function(weekly_age, Scale, shape1, Shift){
  
  dist = (((shape1/Scale)*((weekly_age-Shift)/Scale)^(shape1-1))/((1+((weekly_age-Shift)/Scale)^shape1)^2))
  if(is.na(dist)){
    dist = 0
  }
  return(dist)
}

#9. Takes the start of the weekly age distribution i.e. from 0 to 259 and this will be used for rescaling coverage for timeliness
start_weekly_dist_shift <- function(week, Scale, shape1, Shift){
  
  # R package required to run function
  library(data.table)
  
  #1. CDF Mid week
  pdf_start_week <- list()
  cdf_start_week <- list()
  Start_Week <- list()
  for(i in 0:(week-1)){
    pdf_start_week[i+1] =  burr_shift_pdf(i, Scale, shape1, Shift)
    Start_Week[i+1] = i
    cdf_start_week[i+1] = sum(unlist(pdf_start_week[1:i+1]))
  }
  start_Week_cdf <- as.data.table(cbind(Start_Week, cdf_start_week))
}


#10. Distribution of i.e. cases by year of age
yearly_scenario_func <- function(scenarios_df, dist_df, scenario_type){
  
  years = unlist(unique(scenarios_df$Y))
  diseases = unlist(unique(scenarios_df$Disease))
  pdf = unlist(dist_df$PDF)
  
  cases_weekly <- list()
  under_1s<- list()
  one <- list()
  two <- list()
  three <- list()
  four <- list()
  Total <- list()
  a <- 0
  Y <- list()
  Disease <- list()
  
  for(j in 1:length(years)){
    per_year_df = scenarios_df[Y == years[j],]
    Scenario = per_year_df[[scenario_type]]
    
    #### YEARLY NO. OF SCENARIO'S FOR EACH COHORT ####
    for(i in 1:length(diseases)){
      a = a + 1
      Y[a] = years[j]
      Disease[a] <- as.character(diseases[i])
      
      scenario_weekly = as.numeric(Scenario[i])*pdf         #make a function to retrieve weekly cases 
      under_1s[[a]] = sum(scenario_weekly[1:52])
      one[[a]] = sum(scenario_weekly[53:104]) 
      two[[a]] = sum(scenario_weekly[105:156])
      three[[a]] = sum(scenario_weekly[157:208])
      four[[a]] = sum(scenario_weekly[209:260])
      Total[[a]] = sum(scenario_weekly)
    } 
  }
  
  Yearly_scenarios <- as.data.table(cbind(Disease, Y, Total, under_1s, one, two, three, four))
  df2 <- as.data.frame(lapply(Yearly_scenarios, unlist))
  final_df <- df2[order(df2[, 1]), ]
  return(final_df)
}


#11. Function which multiplies cdf and pdf of total cases/scenarios to get the no. of mid-weekly cases at week x - age distribution
weekly_scenario_func <- function(scenarios_df, dist_df, scenario_type){
  
  years = unlist(unique(scenarios_df$Y))
  diseases = unlist(unique(scenarios_df$Disease))
  pdf = unlist(dist_df$PDF)
  cdf = unlist(dist_df$RESCALED_CDF)
  
  cases_weekly <- list()
  Total <- list()
  a <- 0
  Y <- list()
  Disease <- list()
  
  for(j in 1:length(years)){
    per_year_df = scenarios_df[Y == years[j],]
    Scenario = per_year_df[[scenario_type]]    #takes the column values for a particular scenario i.e. 'Cases'
    
    #### YEARLY NO. OF SCENARIO'S FOR EACH COHORT ####
    for(i in 1:length(diseases)){
      a = a + 1
      Y[a] = years[j]
      Disease[a] <- as.character(diseases[i])
      
      cumulative_weekly_scenario = as.numeric(Scenario[i])*cdf
      scenario_weekly = as.numeric(Scenario[i])*pdf         #make a function to retrieve weekly cases
      Total[[a]] = sum(scenario_weekly)
      cases_weekly[[a]] = rbind(scenario_weekly)
      
    } 
  }
  m = matrix(unlist(cases_weekly), ncol = 260, byrow = TRUE)
  Weekly_scenarios <- as.data.table(cbind(Disease, Y, Total, m))
  df2 <- as.data.frame(lapply(Weekly_scenarios, unlist))
  final_df <- df2[order(df2[, 1]), ]
  colnames(final_df)[4:263] <- as.character(c(1:260))
  
  return(final_df)
}


#12. Dalys before vaccination
Dalys_preVax <- function(NoVax_cases, NoVax_deaths_yearly, ERYOL, DALYs_dt, run.id, random.no, probability, Random.no.list, psa){
  
  NoVax_cases = as.matrix(NoVax_cases)
  NoVax_deaths_yearly = as.data.table(NoVax_deaths_yearly)
  
  year = as.numeric(NoVax_cases[2])
  disease_type = NoVax_cases[1]
  NoVax_deaths = NoVax_deaths_yearly[Y == year & Disease == disease_type]
  
  Exp_RYOL = ERYOL[cohort_year == year]
  DALYs = DALYs_dt[disease == disease_type]
  country_code = Exp_RYOL[,country]
  pathogen_name = DALYs$pathogen
  
  
  if(psa == T){
    # work out DALYs by first calculating YLDs and YLLs for each cohort year
    # YLD = no. of cases by given age x healthy time lost for given disease x duration of infection for given disease
    healthy_time_lost_Mid = DALYs$Mid_healthy_time_lost
    healthy_time_lost_Low = DALYs$Low_healthy_time_lost
    healthy_time_lost_High = DALYs$High_healthy_time_lost
    healthy_time_lost = Beta.Pert.dist(healthy_time_lost_Mid, healthy_time_lost_Low, healthy_time_lost_High, probability)
    
    healthy_time_lost_PSA <- data.frame()
    if(year == 1996){
      healthy_time_lost_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, healthy_time_lost))
      
    }
    
    random.no = random.no + 1
    probability = Random.no.list[random.no]
    duration_of_illness_Mid = DALYs$Mid_duration_of_illness
    duration_of_illness_Low = DALYs$Low_duration_of_illness
    duration_of_illness_High = DALYs$High_duration_of_illness
    duration_of_illness = Beta.Pert.dist(duration_of_illness_Mid, duration_of_illness_Low, duration_of_illness_High, probability)
    
    duration_of_illness_PSA <- data.frame()
    if(year == 1996){
      duration_of_illness_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, duration_of_illness))
      
    }
  }else if(psa == F){
    healthy_time_lost_PSA <- data.frame()
    healthy_time_lost = DALYs$Mid_healthy_time_lost
    
    if(year == 1996){
      healthy_time_lost_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, healthy_time_lost))
      
    }
    random.no = random.no + 1
    probability = Random.no.list[random.no]
    duration_of_illness = DALYs$Mid_duration_of_illness
    duration_of_illness_PSA <- data.frame()
    if(year == 1996){
      duration_of_illness_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, duration_of_illness))
      
    }
  }
  
  YLD_0 = as.numeric(NoVax_cases[4]) * healthy_time_lost * duration_of_illness
  YLD_1 = as.numeric(NoVax_cases[5]) * healthy_time_lost * duration_of_illness
  YLD_2 = as.numeric(NoVax_cases[6]) * healthy_time_lost * duration_of_illness
  YLD_3 = as.numeric(NoVax_cases[7]) * healthy_time_lost * duration_of_illness
  YLD_4 = as.numeric(NoVax_cases[8]) * healthy_time_lost * duration_of_illness
  
  # YLL = no.of deaths for given year x expected remaining years of life by given year and cohort
  YLL_0 = as.numeric(NoVax_deaths[, under_1s]) * Exp_RYOL$X0
  YLL_1 = as.numeric(NoVax_deaths[, one]) * Exp_RYOL$X1
  YLL_2 = as.numeric(NoVax_deaths[, two]) * Exp_RYOL$X2
  YLL_3 = as.numeric(NoVax_deaths[, three]) * Exp_RYOL$X3
  YLL_4 = as.numeric(NoVax_deaths[, four]) * Exp_RYOL$X4
  
  # DALYs = YLD + YLL
  DALY_0 = YLD_0 + YLL_0
  DALY_1 = YLD_1 + YLL_1
  DALY_2 = YLD_2 + YLL_2
  DALY_3 = YLD_3 + YLL_3
  DALY_4 = YLD_4 + YLL_4
  
  dalys_preVax = c(country_code, disease_type, year, DALY_0, DALY_1, DALY_2, DALY_3, DALY_4)
  names(dalys_preVax) <- c("country_code", "disease", "year", "dalys_0", "dalys_1", "dalys_2", "dalys_3", "dalys_4")
  #return(dalys_preVax)
  mydata <- list(dalys_preVax = dalys_preVax, healthy_time_lost_PSA = healthy_time_lost_PSA, duration_of_illness_PSA = duration_of_illness_PSA)
  return(mydata)
  
}


#13. imapact calculation
Vax_impact <- function(pathogen_name, run.id, random.no, probability, Random.no.list, country_code, country, NoVax_Cases,
                       NoVax_deaths_weekly,country_pop,vaccine_start=2000,
                       country_coverage,d1_timeliness, d2_timeliness, d3_timeliness,
                       disease_efficacy, efficacy_dt, ERYOL, DALYs_dt,impact_table,psa, waning){

  library(data.table)
  NoVax_deaths_weekly = as.data.table(NoVax_deaths_weekly)
  year = as.numeric(NoVax_Cases[2])
  disease_type = as.matrix(NoVax_Cases)[1]
  Cases_NoVax = as.numeric(NoVax_Cases[4:length(NoVax_Cases)])
  NoVax_Deaths = NoVax_deaths_weekly[Y == year & Disease == disease_type]
  Deaths_NoVax = as.numeric(NoVax_Deaths[1,4:length(NoVax_Deaths)])

  #use mid-weekly age converted to months to calculate gamma distribution for efficacy waning
  mid_weeks_to_months = ((0.5:259.50)*7)/(365.25/12)

  #getting the efficacy and paramters for the gamma distribution by Dose and Efficacy type (E1 or E2)
  efficacy_type = disease_efficacy[disease == disease_type]$efficacy

  #if(!(efficacy_type %in% impact_table$V1 & year %in% impact_table$V2)){
  if(nrow(impact_table[V1 == efficacy_type & V2 == year]) == 0) {

    efficacy_d1_params= efficacy_dt[pathogen == pathogen_name & Dose == 'D1' & efficacy == efficacy_type]
    efficacy_d2_params= efficacy_dt[pathogen == pathogen_name & Dose == 'D2' & efficacy == efficacy_type]
    efficacy_d3_params= efficacy_dt[pathogen == pathogen_name & Dose == 'D3' & efficacy == efficacy_type]
    
    if(psa == T){
      efficacy_d1_mid = efficacy_d1_params[gamma_params == 'VE', mid]
      efficacy_d1_low = efficacy_d1_params[gamma_params == 'VE', low]
      efficacy_d1_high = efficacy_d1_params[gamma_params == 'VE', high]
      efficacy_d1 = Beta.Pert.dist(efficacy_d1_mid, efficacy_d1_low, efficacy_d1_high, probability)
      
      efficacy_d2_mid = efficacy_d2_params[gamma_params == 'VE', mid]
      efficacy_d2_low = efficacy_d2_params[gamma_params == 'VE', low]
      efficacy_d2_high = efficacy_d2_params[gamma_params == 'VE', high]
      efficacy_d2 = Beta.Pert.dist(efficacy_d2_mid, efficacy_d2_low, efficacy_d2_high, probability)
      
      efficacy_d3_mid = efficacy_d3_params[gamma_params == 'VE', mid]
      efficacy_d3_low = efficacy_d3_params[gamma_params == 'VE', low]
      efficacy_d3_high = efficacy_d3_params[gamma_params == 'VE', high]
      efficacy_d3 = Beta.Pert.dist(efficacy_d3_mid, efficacy_d3_low, efficacy_d3_high, probability)
    }else if(psa == F){
      
      efficacy_d1 = efficacy_d1_params[gamma_params == 'VE', mid]
      efficacy_d2 = efficacy_d2_params[gamma_params == 'VE', mid]
      efficacy_d3 = efficacy_d3_params[gamma_params == 'VE', mid]
      
    }

    if(year == 2000){

      efficacy_PSA <- data.frame()
      efficacy_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, efficacy_type, efficacy_d1, efficacy_d2, efficacy_d3))

    }


    random.no = random.no + 1
    probability = Random.no.list[random.no]
    #------------------------------------------------------------------------------------------------------
    
    if(psa == T){
      
      Alpha_d1_mid = efficacy_d1_params[gamma_params == 'alpha', mid]
      Alpha_d1_low = efficacy_d1_params[gamma_params == 'alpha', low]
      Alpha_d1_high = efficacy_d1_params[gamma_params == 'alpha', high]
      Alpha_d1 = Beta.Pert.dist(Alpha_d1_mid, Alpha_d1_low, Alpha_d1_high, probability)
      
      Mean_d1_mid = efficacy_d1_params[gamma_params == 'Mean', mid]
      Mean_d1_low = efficacy_d1_params[gamma_params == 'Mean', low]
      Mean_d1_high = efficacy_d1_params[gamma_params == 'Mean', high]
      Mean_d1 = Beta.Pert.dist(Mean_d1_mid, Mean_d1_low, Mean_d1_high, probability)
      
      #------------------------------------------------------------------------------------------------------
      
      Alpha_d2_mid = efficacy_d2_params[gamma_params == 'alpha', mid]
      Alpha_d2_low = efficacy_d2_params[gamma_params == 'alpha', low]
      Alpha_d2_high = efficacy_d2_params[gamma_params == 'alpha', high]
      Alpha_d2 = Beta.Pert.dist(Alpha_d2_mid, Alpha_d2_low, Alpha_d2_high, probability)
      
      Mean_d2_mid = efficacy_d2_params[gamma_params == 'Mean', mid]
      Mean_d2_low = efficacy_d2_params[gamma_params == 'Mean', low]
      Mean_d2_high = efficacy_d2_params[gamma_params == 'Mean', high]
      Mean_d2 = Beta.Pert.dist(Mean_d2_mid, Mean_d2_low, Mean_d2_high, probability)
      
      #------------------------------------------------------------------------------------------------------
      
      Alpha_d3_mid = efficacy_d3_params[gamma_params == 'alpha', mid]
      Alpha_d3_low = efficacy_d3_params[gamma_params == 'alpha', low]
      Alpha_d3_high = efficacy_d3_params[gamma_params == 'alpha', high]
      Alpha_d3 = Beta.Pert.dist(Alpha_d3_mid, Alpha_d3_low, Alpha_d3_high, probability)
      
      Mean_d3_mid = efficacy_d3_params[gamma_params == 'Mean', mid]
      Mean_d3_low = efficacy_d3_params[gamma_params == 'Mean', low]
      Mean_d3_high = efficacy_d3_params[gamma_params == 'Mean', high]
      Mean_d3 = Beta.Pert.dist(Mean_d3_mid, Mean_d3_low, Mean_d3_high, probability)
      
    }else if(psa == F){
      
      if(waning == T){
        Alpha_d1 = efficacy_d1_params[gamma_params == 'alpha', mid]
        Mean_d1 = efficacy_d1_params[gamma_params == 'Mean', mid]
        
        Alpha_d2 = efficacy_d2_params[gamma_params == 'alpha', mid]
        Mean_d2 = efficacy_d2_params[gamma_params == 'Mean', mid]
        
        Alpha_d3 = efficacy_d3_params[gamma_params == 'alpha', mid]
        Mean_d3 = efficacy_d3_params[gamma_params == 'Mean', mid]
      }else if(waning == F){
        Alpha_d1 = 100
        Mean_d1  = 1000
        
        Alpha_d2 = 100
        Mean_d2  = 1000
        
        Alpha_d3 = 100
        Mean_d3  = 1000
      }
      
    }

    Beta_d1 = Mean_d1/Alpha_d1
    gamma_distribution_d1 = sapply(mid_weeks_to_months, function(x) 1- pgamma(x, shape = Alpha_d1, scale = Beta_d1))
    rescale_VE_2wks_d1 = 1 - pgamma(0.5, shape = Alpha_d1, scale = Beta_d1)

    Beta_d2 = Mean_d2/Alpha_d2
    gamma_distribution_d2 = sapply(mid_weeks_to_months, function(x) 1- pgamma(x, shape = Alpha_d2, scale = Beta_d2))
    rescale_VE_2wks_d2 = 1 - pgamma(0.5, shape = Alpha_d2, scale = Beta_d2)

    Beta_d3 = Mean_d3/Alpha_d3
    gamma_distribution_d3 = sapply(mid_weeks_to_months, function(x) 1- pgamma(x, shape = Alpha_d3, scale = Beta_d3))
    rescale_VE_2wks_d3 = 1 - pgamma(0.5, shape = Alpha_d3, scale = Beta_d3)

    #------------------------------------------------------------------------------------------------------

    if(year == 2000){

      efficacy_PSA2 <- data.frame()
      efficacy_PSA2 <- rbind(cbind(country, run, random.no, probability, pathogen_name, efficacy_type, Mean_d1, Mean_d2, Mean_d3, Alpha_d1, Alpha_d2, Alpha_d3, Beta_d1, Beta_d2, Beta_d3))

    }

    #2. Get the coverage data for each year and then the efficacy for each dose, #for timeliness use exact coverage
    path_name = pathogen_name
    if(path_name == "Rota"){
      path_name = "RV"}

    if(as.numeric(NoVax_Cases[2]) >= vaccine_start){
      cov_year = country_coverage[YEAR == as.numeric(NoVax_Cases[2])]
      d1 = paste(path_name,"_DTP1",sep="")
      d2 = paste(path_name,"_DTP2",sep="")
      d3 = paste(path_name,"_DTP3",sep="")

      cov_d1 = as.numeric(cov_year[, ..d1])
      cov_d2 = as.numeric(cov_year[, ..d2])
      cov_d3 = as.numeric(cov_year[, ..d3])
    } else{cov_year = as.numeric(NoVax_Cases[2])
    cov_d1 = 0
    cov_d2 = 0
    cov_d3 = 0}


    if(cov_d1 == 0 & cov_d2 == 0 & cov_d3 == 0){
      total_impact = rep(0, 260)

      with_vax_cases = Cases_NoVax
      total_cases = sum(with_vax_cases)
      Cases_impact = c(disease_type, year, total_cases, with_vax_cases)


      with_vax_deaths = Deaths_NoVax
      total_deaths = sum(with_vax_deaths)
      Deaths_impact = c(disease_type, year, total_deaths, with_vax_deaths)

      append_row = c(efficacy_type, year, total_impact)
      impact_table = rbind(impact_table, t(append_row))
      #list2env(impact_table, envir = .GlobalEnv)
      assign("impact_table", impact_table, envir = .GlobalEnv)

    }else{

      Mid_Week = seq(0.5, 259.5, 1)
      rescaled_d1_cov = rescale_to_coverage(unlist(d1_timeliness$cdf_start_week), cov_d1, Mid_Week)
      rescaled_d2_cov = rescale_to_coverage(unlist(d2_timeliness$cdf_start_week), cov_d2, Mid_Week)
      rescaled_d3_cov = rescale_to_coverage(unlist(d3_timeliness$cdf_start_week), cov_d3, Mid_Week)
      
      
      #3. Applying waning gamma distribution to efficacy -- HAVE option to turn waning off here...
      d1_efficacy_waning = efficacy_d1/rescale_VE_2wks_d1 * gamma_distribution_d1 
      d2_efficacy_waning = efficacy_d2/rescale_VE_2wks_d2 * gamma_distribution_d2
      d3_efficacy_waning = efficacy_d3/rescale_VE_2wks_d3 * gamma_distribution_d3
      
      
      #4. Get the cumulative cov of each dose and do a Validity check to make sure preceeding coverage of a higher dose is not greater than the previous coverage of previous dose
      cumulative_cov_by_age_mat = as.matrix(rbind(unlist(rescaled_d1_cov$RESCALED_CDF), unlist(rescaled_d2_cov$RESCALED_CDF), unlist(rescaled_d3_cov$RESCALED_CDF)))
      cumulative_cov_by_age = cumulative_covC(cumulative_cov_by_age_mat)
      
      #5. Calculate the incremental coverage by subtraction between doses
      
      incremental_cov_by_age = incremental_covC(cumulative_cov_by_age)

      if(is.element('TRUE', incremental_cov_by_age[1,]<0)){
        print("Stop! Negative incremental coverage in Dose 1")
      }else if(is.element('TRUE', incremental_cov_by_age[3,]<0)){
        print("Stop! Negative incremental coverage in Dose 2")
      }else if(is.element('TRUE', incremental_cov_by_age[3,]<0)){
        print("Stop! Negative incremental coverage in Dose 3")
      }

      #6. Calculate the incremental coverage of each week and time essentially putting in a diagonal matrix
      #d1_weekly_cov in the incremental coverage by week and is organised diagnonally on the matrix
      d1_weekly_cov = diag(incremental_cov_by_age[1,])
      d1_weekly_cov[d1_weekly_cov == 0] = NA

      d2_weekly_cov = diag(incremental_cov_by_age[2,])
      d2_weekly_cov[d2_weekly_cov == 0] = NA

      d3_weekly_cov = as.matrix(diag(incremental_cov_by_age[3,]))
      d3_weekly_cov[d3_weekly_cov == 0] = NA
      d3_weekly_cov[cbind(c(1,2,3),c(1,2,3))] = incremental_cov_by_age[3,1:3]


      #7. Drag the incremental coverage for each week across time
      incremental_cov_d1 = (t(sapply(1:nrow(d1_weekly_cov), function(x){ c(rep(NA, x-1), rep(sum(d1_weekly_cov[x,], na.rm=T),ncol(d1_weekly_cov)-(x-1))) })))
      incremental_cov_d2 = (t(sapply(1:nrow(d2_weekly_cov), function(x){ c(rep(NA, x-1), rep(sum(d2_weekly_cov[x,], na.rm=T),ncol(d2_weekly_cov)-(x-1))) })))
      incremental_cov_d3 = (t(sapply(1:nrow(d3_weekly_cov), function(x){ c(rep(NA, x-1), rep(sum(d3_weekly_cov[x,], na.rm=T),ncol(d3_weekly_cov)-(x-1))) })))


      #8. Calculate the leakage of incremental coverage between doses
      D1_incremental_cov_leaked <- incremental_cov_leakage(incremental_cov_d1, incremental_cov_d2)
      D2_incremental_cov_leaked <- incremental_cov_leakage(incremental_cov_d2, incremental_cov_d3)


      #9. Sum the leaked incremental coverage over weeks and time Additional table to show that the cumulative coverage entered originally can be retrieved by
      D1_incremental_cov_summed = colSums(D1_incremental_cov_leaked, na.rm = T)
      D2_incremental_cov_summed = colSums(D2_incremental_cov_leaked, na.rm = T)
      D3_incremental_cov_summed = colSums(incremental_cov_d3, na.rm = T)

      #10. Additional table to show that the sum of all leaked incremental coverages equates to the cumulative cov for Dose 1
      incremental_cov_leaked_table = rbind(D1_incremental_cov_summed, D2_incremental_cov_summed, D3_incremental_cov_summed)
      incremental_cov_for_all_doses_summed = colSums(incremental_cov_leaked_table)

      #11. Expected impact - incremental_cov[t] * d1_efficacy_waning[t]
      impact_tracked <- impact_time_week(D1_incremental_cov_leaked, D2_incremental_cov_leaked, incremental_cov_d3, d1_efficacy_waning, d2_efficacy_waning, d3_efficacy_waning)
      total_impact = colSums(impact_tracked, na.rm = T)

      append_row = c(efficacy_type, year, total_impact)
      impact_table = rbind(impact_table, t(append_row))
      #list2env(impact_table, envir = .GlobalEnv)
      assign("impact_table", impact_table, envir = .GlobalEnv)

      with_vax_cases <- c()
      for(i in 1:length(total_impact)){
        with_vax_cases[i] = Cases_NoVax[i]*(1 - total_impact[i])
      }

      total_cases = sum(with_vax_cases)
      Cases_impact = c(disease_type, year, total_cases, with_vax_cases)

      with_vax_deaths <- c()
      for(i in 1:length(total_impact)){
        with_vax_deaths[i] = Deaths_NoVax[i]*(1 - total_impact[i])
      }

      total_deaths = sum(with_vax_deaths)
      Deaths_impact = c(disease_type, year, total_deaths, with_vax_deaths)
    }

  }else{

    total_impact = as.numeric(impact_table[V1 == efficacy_type &  V2 == year, 3:262])

    with_vax_cases <- c()
    for(i in 1:length(total_impact)){
      with_vax_cases[i] = Cases_NoVax[i]*(1 - total_impact[i])
    }

    total_cases = sum(with_vax_cases)
    Cases_impact = c(disease_type, year, total_cases, with_vax_cases)

    with_vax_deaths <- c()
    for(i in 1:length(total_impact)){
      with_vax_deaths[i] = Deaths_NoVax[i]*(1 - total_impact[i])
    }

    total_deaths = sum(with_vax_deaths)
    Deaths_impact = c(disease_type, year, total_deaths, with_vax_deaths)

  }



  cohort_pop = country_pop[Year == year]
  Exp_RYOL = ERYOL[cohort_year == year]
  DALYs = DALYs_dt[pathogen == pathogen_name & disease == disease_type]

  zero_cases = sum(with_vax_cases[1:52])
  one_cases = sum(with_vax_cases[53:104])
  two_cases = sum(with_vax_cases[105:156])
  three_cases = sum(with_vax_cases[157:208])
  four_cases = sum(with_vax_cases[209:260])
  cases_by_age = c(zero_cases, one_cases, two_cases, three_cases, four_cases)

  zero_deaths = sum(with_vax_deaths[1:52])
  one_deaths = sum(with_vax_deaths[53:104])
  two_deaths = sum(with_vax_deaths[105:156])
  three_deaths = sum(with_vax_deaths[157:208])
  four_deaths = sum(with_vax_deaths[209:260])
  deaths_by_age = c(zero_deaths, one_deaths, two_deaths, three_deaths, four_deaths)

  # work out DALYs by first calculating YLDs and YLLs for each cohort year
  # YLD = no. of cases by given age x healthy time lost for given disease x duration of infection for given disease

  random.no = random.no + 1
  probability = Random.no.list[random.no]
  
  if(psa == T){
    
    healthy_time_lost_Mid = DALYs$Mid_healthy_time_lost
    healthy_time_lost_Low = DALYs$Low_healthy_time_lost
    healthy_time_lost_High = DALYs$High_healthy_time_lost
    healthy_time_lost = Beta.Pert.dist(healthy_time_lost_Mid, healthy_time_lost_Low, healthy_time_lost_High, probability)
    
    healthy_time_lost_PSA <- data.frame()
    if(year == 2000){
      healthy_time_lost_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, healthy_time_lost))
      
    }
    
    random.no = random.no + 1
    probability = Random.no.list[random.no]
    duration_of_illness_Mid = DALYs$Mid_duration_of_illness
    duration_of_illness_Low = DALYs$Low_duration_of_illness
    duration_of_illness_High = DALYs$High_duration_of_illness
    duration_of_illness = Beta.Pert.dist(duration_of_illness_Mid, duration_of_illness_Low, duration_of_illness_High, probability)
    
    duration_of_illness_PSA <- data.frame()
    if(year == 2000){
      duration_of_illness_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, duration_of_illness))
      
    }
    
  }else if(psa == F){
    
    healthy_time_lost = DALYs$Mid_healthy_time_lost
    healthy_time_lost_PSA <- data.frame()
    if(year == 2000){
      healthy_time_lost_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, healthy_time_lost))
      
    }
    
    random.no = random.no + 1
    probability = Random.no.list[random.no]
    duration_of_illness = DALYs$Mid_duration_of_illness
    duration_of_illness_PSA <- data.frame()
    if(year == 2000){
      duration_of_illness_PSA <- rbind(cbind(country, run, random.no, probability, pathogen_name, disease_type, duration_of_illness))
      
    }
    
  }

  YLD_0 = zero_cases *   healthy_time_lost * duration_of_illness
  YLD_1 = one_cases *    healthy_time_lost * duration_of_illness
  YLD_2 = two_cases *    healthy_time_lost * duration_of_illness
  YLD_3 = three_cases *  healthy_time_lost * duration_of_illness
  YLD_4 = four_cases *   healthy_time_lost * duration_of_illness

  # YLL = no.of deaths for given year x expected remaining years of life by given year and cohort
  YLL_0 = zero_deaths * Exp_RYOL$X0
  YLL_1 = one_deaths * Exp_RYOL$X1
  YLL_2 = two_deaths * Exp_RYOL$X2
  YLL_3 = three_deaths * Exp_RYOL$X3
  YLL_4 = four_deaths * Exp_RYOL$X4

  # DALYs = YLD + YLL
  DALY_0 = YLD_0 + YLL_0
  DALY_1 = YLD_1 + YLL_1
  DALY_2 = YLD_2 + YLL_2
  DALY_3 = YLD_3 + YLL_3
  DALY_4 = YLD_4 + YLL_4
  dalys_by_age = c(DALY_0, DALY_1, DALY_2, DALY_3, DALY_4)

  combined_output <- c(disease_type, year, country_code, as.character(cohort_pop$Country), c(unlist(cohort_pop[,X0:X4])), cases_by_age, deaths_by_age, dalys_by_age)
  names(combined_output) <- c("disease", "year", "country", "country_name", "cohort_0", "cohort_1", "cohort_2", "cohort_3", "cohort_4",
                              "cases_0", "cases_1", "cases_2", "cases_3", "cases_4", "deaths_0", "deaths_1", "deaths_2", "deaths_3", "deaths_4",
                              "dalys_0", "dalys_1", "dalys_2", "dalys_3", "dalys_4")

  mydata <- list(Cases = Cases_impact, Deaths = Deaths_impact, combined_output = combined_output, efficacy_PSA_output = efficacy_PSA, efficacyPSA2 = efficacy_PSA2,
                 healthy_time_lost_PSA = healthy_time_lost_PSA, duration_of_illness_PSA = duration_of_illness_PSA)

  return(mydata)


}
