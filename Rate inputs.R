#------------------------------------------------------------------
# Title: Rate inputs
# Author: Eleanor Hayes-Larson
# Co-author/editor: Crystal Shaw
# Purpose: This program loads in data to be used in calibrating the 
#   multistate model.
#------------------------------------------------------------------

#------------------------------------------------------------------
#---- Defining mortality rates ----
#------------------------------------------------------------------
#Make overall mortality rate a function of time
#---- Load survival data from lifetables ----
Lifetables <- read.csv(file = here("Calibration data","Lifetables.csv"),
                       header = T, sep = ",", stringsAsFactors = F) %>%
  dplyr::select(c("Age", "Allrace_all", "Allrace_M", "Allrace_F")) 

#---- Calculate conditional mortality rate and cumulative survival ----
#Mort rate is # deaths/person-time 
#(5 years for survivors and 2.5 for deaths in interval)

#Create empty columns for mortality rates overall and by sex/gender
genders <- c("all", "M", "F")
Lifetables %<>% 
  cbind(., as.data.frame(matrix(nrow = nrow(Lifetables), ncol = 3)) %>% 
          set_colnames(paste0("Mortality_obs_", genders))) %>% 
  cbind(., as.data.frame(matrix(nrow = nrow(Lifetables), ncol = 3)) %>% 
          set_colnames(paste0("Cum_surv_obs_", genders)))

#Fill in the table by gender
for(gender in genders){
  #Mortality rate
  Allrace_column <- paste0("Allrace_", gender) #Define the column name
  
  Lifetables[, paste0("Mortality_obs_", gender)] <- 
    (lag(Lifetables[, Allrace_column]) - Lifetables[, Allrace_column])/
    (5*Lifetables[, Allrace_column] + 
       2.5*(lag(Lifetables[, Allrace_column]) - Lifetables[, Allrace_column]))
  
  #Cumulative survival
  Lifetables[, paste0("Cum_surv_obs_", gender)] <- 
    (Lifetables[, Allrace_column]/Lifetables[, Allrace_column][1])
}

Lifetables$time <- seq(from = 0, to = 40, by = 5)
Lifetables

#------------------------------------------------------------------
#---- Defining cancer incidence rates ---- 
#------------------------------------------------------------------
#make cancer incidence rate a function of time 
#(study time = age, starting at age 65) 
ages <- c("65-69", "70-74", "75-79", "80-84", "85+")
types <- c("all", "lung", "breast", "prostate")

SEER_inc_obs <- 
  data.frame(Age = ages, time = seq(from = 5, to = 25, by = 5)) %>%
  cbind(., as.data.frame(matrix(nrow = length(ages), 
                                ncol = length(types))) %>% 
          set_colnames(paste0("Rate_", types)))
                    
for (type in types){
  SEER_inc <- read.csv(here("Calibration data", 
                            paste0("SEER ", type, " cancer incidence.csv")),
                       skip = 3, stringsAsFactors = F, header = T)
  
  if(type == "all" | type == "lung"){
    sex = "Both Sexes"
  } else if(type == "breast"){
    sex = "Female"
  } else {sex = "Male"}
  
  #Pull rates from SEER, and convert from rate per 100k.
    if (type != "lung"){
      temp_data <- SEER_inc %>% 
        filter(Age %in% ages, Sex == sex, 
               Race.Ethnicity == "All Races (includes Hispanic)",  
               Rate.Type == "Delay-adjusted Rates") %>% 
        dplyr::select("Rate.per.100.000") %>% 
        mutate_all(~as.numeric(.)/100000) 
      SEER_inc_obs[, paste0("Rate_", type)] <- temp_data
    } else {
      #Lung cancer SEER data doesn't have a label for observed vs. 
      #delay-adjusted, but larger rates (rows 6-10) are delay-adjusted, 
      #so I take these.
      temp_data <- SEER_inc %>% 
        filter(Age %in% ages, Sex == sex,  
               Race.Ethnicity == "All Races (includes Hispanic)") %>% 
        dplyr::select("Rate.per.100.000") %>% mutate_all(~as.numeric(.)/100000)
      SEER_inc_obs[, paste0("Rate_", type)] <- temp_data[6:10, ]
    } 
}

SEER_inc_obs

#------------------------------------------------------------------
# ---- Defining cancer relative survival ----
#------------------------------------------------------------------
#make cancer relative survival a function of time 
#(study time = age, starting at age 65) 

ages2 <- c("65-74", "75+")
SEER_rel.surv <- data.frame(Age = ages2) %>% 
  cbind(., as.data.frame(matrix(nrow = length(ages2), 
                                ncol = length(types))) %>% 
          set_colnames(paste0("Rel.surv_", types)))
                          

for (type in types){
  SEER_surv <- read.csv(here("Calibration data", 
                             paste0("SEER ", type, " cancer survival.csv")),
                        skip = 3, stringsAsFactors = F, header = T)
  
  if(type == "all" | type == "lung"){
    sex = "Both Sexes"
  } else if(type == "breast"){
    sex = "Female"
  } else {sex = "Male"}
  
  Rel.surv_column <- paste0("Rel.surv_", type) #Define the column name
  
  temp_data <- SEER_surv %>% 
    filter(Age %in% paste0("Ages ", ages2), Sex == sex, 
           Race.Ethnicity == "All Races (includes Hispanic)",  
           Stage.at.Diagnosis == "All Stages") %>% 
    dplyr::select("Relative.Survival.Rate....") %>% 
    mutate_all(~as.numeric(.)/100) #convert to decimal
    
  SEER_rel.surv[, Rel.surv_column] <- temp_data
}

SEER_rel.surv

#------------------------------------------------------------------
# ---- Defining dementia incidence rates ----
#------------------------------------------------------------------
#make dementia incidence rate a function of time 
#(study time = age, starting at age 65)
#Load ACT study data--use all-cause dementia, not just AD.
AD_inc_obs <- read.csv(file = here("Calibration data","AD_inc_Tom_etal.csv"),
                       header = T, sep = ",", stringsAsFactors = F) 

#Convert from rate/1000
rate_columns <- str_detect(colnames(AD_inc_obs), "Rate")
AD_inc_obs[, rate_columns] <- AD_inc_obs[, rate_columns]/1000

AD_inc_obs$time <- seq(from = 0, to = 40, by = 5)

AD_inc_obs

#------------------------------------------------------------------
# ---- Defining dementia mortality hazard ratios ----
#------------------------------------------------------------------
#Load Mayeda 2017 study data--use results for whites
AD_HR_obs <- read.csv(file = here("Calibration data", 
                                  "Dem_mort_Mayeda_etal.csv"),
                      header = T, sep = ",", stringsAsFactors = F)
AD_HR_obs


#------------------------------------------------------------------
# ---- Dropping unneeded objects ----
#------------------------------------------------------------------

rm(list = setdiff(ls(), 
                  c("Lifetables", "SEER_inc_obs", "SEER_rel.surv", "AD_inc_obs",
                    "AD_HR_obs",lsf.str())))
