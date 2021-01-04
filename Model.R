#------------------------------------------------------------------
# Title:Compartment model function
# Author: Eleanor Hayes-Larson
# Co-author/Editor: Crystal Shaw
#------------------------------------------------------------------

#Input parameters of simulation
#immortal #This is a switch to make cohort immortal (when immortal == 1)
#diff.cancer.mort #This is a switch to make mortality higher in those with cancer
#diff.dem.mort #This is a switch to make mortality higher in those with dementia
#S1.dem.rateratio #Annual rate ratio for protective effect of S charactistic on dementia
#S1.mort.rateratio #Annual rate ratio for protective effect of S characteristic on mortality in cancer+
#times  #Time steps for model in years. Request long time interval to ensure steady state (100% dead)

#---- Function ----
compartmodel <- function(immortal, diff.cancer.mort, diff.dem.mort,
                         S1.dem.rateratio, S1.mort.rateratio, C1D1mortmultiplier, 
                         Sprev, 
                         Incr_contact, Decr_demdx, Demdxwaittime,
                         times, 
                         type){
  
  if(type == "all" | type == "lung"){
    sex = "all"
  } else if(type == "breast"){
    sex = "F"
  } else {
    sex = "M"
  }
  
  
  Lifetab <- Lifetables[, c("Age", paste0("Mortality_obs_", sex))] 
  Surv_formod<-Lifetab[2:9,]
  Surv_formod$Mortality_obs<-Surv_formod[,2]
  Surv_formod$Age<-Surv_formod$Age-67.5
  Mortality<-rationalfit(Surv_formod$Age, Surv_formod$Mortality_obs, d1=2,d2=1)
  
  SEER_rel <- SEER_rel.surv[, c("Age", paste0("Rel.surv_", type))] 
  SEER_relsurv_mod<-data.frame(Age=seq(2.5,37.5,5),
                               Rel.surv=c(rep(SEER_rel[SEER_rel$Age=="65-74",2],2),
                                          rep(SEER_rel[SEER_rel$Age=="75+",2],6)))
  Rel_surv<-rationalfit(SEER_relsurv_mod$Age, SEER_relsurv_mod$Rel.surv, d1=3,d2=0)
  
  ADHR_formod<-AD_HR_obs[,c("Age","Dementia.HR")]
  ADHR_formod$Age<-c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5)
  AD_HR_obs<-rationalfit(ADHR_formod$Age, ADHR_formod$Dementia.HR, d1=1,d2=1)
 
  
  SEER_in_formod<- SEER_inc_obs[, c("time",paste0("Rate_", type))]
  SEER_in_formod$Rate<-SEER_in_formod[,2]
  SEER_in_formod$time<-(SEER_in_formod$time-2.5)
  Cancer_inc<-rationalfit(SEER_in_formod$time, SEER_in_formod$Rate, d1=2,d2=0)

  AD_in_formod <- AD_inc_obs[2:9, c("Age",paste0("Rate_", sex))]
  AD_in_formod$Rate<-AD_in_formod[,2]
  AD_in_formod$Age<-(AD_in_formod$Age-67.5)
  Dem_inc<-rationalfit(AD_in_formod$Age, AD_in_formod$Rate, d1=3,d2=1)
  
  
  Demdx_inc <- rep(1/Demdxwaittime,9)
  
  calibration_data <- 
    read_csv(here("Code", "Multistate model", 
                  "calibrated_pars",
                  paste0(type, "_calibration.csv"))) %>% as.matrix()
  
  # Mult<-data.frame(Age=seq(1,35),
  #                  Multi=c(rep(0.43,20),rep(0.70,15)))
  # 
  # Multfactor<-rationalfit(Mult$Age, Mult$Multi, d1=1,d2=0)

  
  #------------------------------------------------------------------
  # # Part 1: Defining sources and calibration parameters for each cancer type
  #   #These are calibrations that I did manually 
  #   #They can be checked using the "Calibration checks" Rmd file.
  #------------------------------------------------------------------
  #   
  # All calibrated parameters are now in their own .csv files and loaded above
  #                                                          
  #------------------------------------------------------------------
  # Part 2: Defining input rates based on time (i.e. age)
  # In this section, I define several functions for rates corresponding to 
  # different states (e.g. C0D1S0). I start with a base rate and apply 
  # rate ratios to get rates for other groups.   
  #------------------------------------------------------------------
  #
  #   These are now all defined in mortality_rates.R script
  #
  #------------------------------------------------------------------
  # Part 4: Run simulations
  #------------------------------------------------------------------
  # lsoda() is a function which numerically solves the differential equation
  # inputs: initial states, times where you want to know the solution, 
  #         the model, and the parameters of the model
  # parms is the list of parameters that other functions in diff.eq.model need:
  ## Cancer_inc, S1.dem.rateratio, Sprev, calibration_data, immortal, 
  ## Mortality, diff.cancer.mort, diff.dem.mort, AD_HR_obs, 
  ## Rel_surv, S1.mort.rateratio
  
  # lsoda automatically determines the solving method and the time step to get 
  # the solution to within some tolerance at the specified times
  
  
  #RUN SIMULATION
  #This runs the simulation and generates the state distribution dataframe
  #Script for differential equation model
  #source(here("Code", "Multistate model", "diff_eq_model.R"))
  
  parms = list("Cancer_inc" = Cancer_inc, "S1.dem.rateratio" = S1.dem.rateratio, 
               "C1D1mortmultiplier" = C1D1mortmultiplier, #Multfactor,
               "Sprev" = Sprev, "calibration_data" = calibration_data, 
               "immortal" = immortal, "Mortality" = Mortality, 
               "diff.cancer.mort" = diff.cancer.mort, 
               "diff.dem.mort" = diff.dem.mort, "AD_HR_obs" = AD_HR_obs, 
               "Rel_surv" = Rel_surv, "S1.mort.rateratio" = S1.mort.rateratio, 
               "Dem_inc" = Dem_inc, "Demdx_inc" = Demdx_inc,
               "Incr_contact" = Incr_contact, "Decr_demdx" = Decr_demdx, 
               "Demdxwaittime"=Demdxwaittime)
  
  out.data <- ode(gen.init(Sprev), times, diff.eq.model, parms, method="ode45")
  out.data %<>% as.data.frame()
  
  #This saves the rate results
  rate.results <- matrix(NA, ncol = length(times), nrow = 43)
  
  for (i in 1:length(times)){
    rate.results[, i] <- 
      diff.eq.model(time = round(times[i], 2), 
                    out.data[round(out.data$time, 2) == round(times[i], 2), ], 
                    parms)[[1]]
  }
  
  rate.results %<>% t() %>% cbind(times, .) %>% 
    set_colnames(c("time", 
                   "dC0D0S0", "dC0D1S0", "dC1D0S0", "dC1D1S0",
                   "dC0D0S1", "dC0D1S1", "dC1D0S1", "dC1D1S1", 
                   "dC0D2S0", "dC0D2S1", "dC1D2S0", "dC1D2S1", 
                   "dDEAD", 
                   "dPT_C0D0S0",
                   "dPT_C0D0S1",
                   "dPT_C1D0S0",
                   "dPT_C1D0S1",
                   "dPT_C0D1S0",
                   "dPT_C0D1S1",
                   "dPT_C1D1S0",
                   "dPT_C1D1S1",
                   "dPT_C0D2S0",
                   "dPT_C0D2S1",
                   "dPT_C1D2S0",
                   "dPT_C1D2S1", 
                   "dDEMNOC", "dDEMWC", "dDEMdxNOC", "dDEMdxWC", 
                   "dDEMENTIA", "dDEMENTIAdx", 
                   "dCANCER", "dDEADC0D0S0", "dDEADC1D0S0", 
                   "dDEADC0D1S0", "dDEADC1D1S0", "dDEADC0D0S1", 
                   "dDEADC1D0S1", "dDEADC0D1S1", "dDEADC1D1S1", 
                   "dCBEFOREDEM", "dDEMNOSNOC", "dDEMWITHSNOC"))
  #This order has to match as above. 
  
  #Outputs for full compartmodel function 
  # out.data is a matrix with each row for the # of people in each state at 
  #input times rate.results is a matrix of rates for each states at each 
  #input time.
  
  return(list(out.data = out.data, rate.results = rate.results))
}

# #---- Part 6: Testing Code ----
# immortal = 0
# diff.cancer.mort = 1
# diff.dem.mort = 1
# S1.dem.rateratio = 0.5
# S1.mort.rateratio = 0.5
# Sprev = 0.5
# times = c(seq(from = 0, to = 40, by = 1/12), 10000)
# type = "all"
# 
# test <- compartmodel(0, 1, 1,
#                      #bin_cuts, IRRs, prevs,
#                      0.5, 0.5, 0.5,
#                      Incr_contact=1, Decr_demdx=1,
#                      c(seq(from = 0, to = 40, by = 1/12), 10000),
#                      "all")
# 

# microbenchmark::microbenchmark(compartmodel(0, 1, 1, bin_cuts, IRRs, prevs,
#                               0.5, 0.5, 0.5,
#                               c(seq(from = 0, to = 40, by = 1/12), 10000),
#                               "all"), times = 10)
