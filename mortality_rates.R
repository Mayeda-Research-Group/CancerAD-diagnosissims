#This script has all of the functions related to the mortality rates

# Surv_formod<-Lifetables[,c("Age","Mortality_obs_all")]
# Surv_formod$Age<-Surv_formod$Age-67.5
# Surv_formod
# test3<-rationalfit(Surv_formod$Age, Surv_formod$Mortality_obs_all, d1=2,d2=4                    )
# p1 <- test3$p1; p2 <- test3$p2
# ys <- polyval(p1,seq(0,40,1/12)) / polyval(p2,seq(0,40,1/12))
# plot(Surv_formod$Age, Surv_formod$Mortality_obs_all, type="l", col = "blue", xlim=c(0,40), ylim=c(0,0.5))
# points(seq(0,40,1/12), Re(ys), col="red")
# points(Surv_formod$Age, Surv_formod$Mortality_obs_all, col="blue", pch=15)

# ADHR_formod<-AD_HR_obs[,c("Age","Dementia.HR")]
# ADHR_formod$Age<-c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5)
# ADHR_formod
# ADHRmod<-rationalfit(ADHR_formod$Age, ADHR_formod$Dementia.HR, d1=2,d2=0                    )
# p1_ADHR <- ADHRmod$p1; p2_ADHR <- ADHRmod$p2
#ys<-(polyval(p1_ADHR,time) / polyval(p2_ADHR,time))
# plot(ADHR_formod$Age, ADHR_formod$Dementia.HR, type="l", col = "blue", xlim=c(0,40), ylim=c(1,10))
# points(seq(0,40,1/12), Re(ys), col="red")
# points(ADHR_formod$Age, ADHR_formod$Dementia.HR, col="blue", pch=15)


# Define reference mortality rate for non-cancer, non-dementia as a 
# function of time using lifetables data for whites
mortality.rate.C0D0 <- function(immortal, S1.dem.rateratio, Sprev, 
                                S1.mort.rateratio, time, calibration_data, 
                                Mortality){
  if (immortal == 1){
    return(0)
  }
  
  Mort_cal <- 
    calibration_data[, paste0("Mort_cal_IRR", S1.dem.rateratio, "_prev", Sprev, 
                              "_mortRR", S1.mort.rateratio)]
  max_index <- length(Mort_cal)
  
  if(time %% 5 == 0){
    index = time/5 + 1
  } else{
    #Rounds time up to the next multiple of 5 then divides by 5 for indexing
    index = ceiling(time/5)
  }
  if(index > max_index){index = max_index}
  # 
  
  Mort_rate <- Re(polyval(Mortality$p1,time) / polyval(Mortality$p2,time))

  return(Re(Mort_rate)*Mort_cal[index])
  #return(Mortality[index+1]*Mort_cal[index])
}

#This function defines mortality IRR for cancer using relative survival
# mort.cancer.irr <- function(time, immortal, S1.dem.rateratio, Sprev, 
#                             S1.mort.rateratio, calibration_data, Mortality, 
#                             rel.surv.5yrs){
#   C0D0_rate <- mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
#                                    S1.mort.rateratio, time, calibration_data, 
#                                    Mortality)
#   
#   
#   mort.cancer.rate = C0D0_rate - 1/5*log(rel.surv.5yrs)
#   
#   return(mort.cancer.rate/C0D0_rate)
# }

#This function defines mortality rates, applying IRR to base rate
mortality.rate.C1D0 <- function(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                Sprev, S1.mort.rateratio, time, 
                                calibration_data, Mortality, Rel_surv){
  if (immortal == 1){return(0)}
  if (diff.cancer.mort == 0){
    return(mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
                               S1.mort.rateratio, time, calibration_data, 
                               Mortality))
  } else{
    
    #Effect of cancer on mortality
    rel.surv.5yrs<- Re(polyval(Rel_surv$p1,time) / polyval(Rel_surv$p2,time))
    ratediff<- as.numeric(- 1/5*log(rel.surv.5yrs))
    ratediff<-ifelse(ratediff<0,0,ratediff)
    
      return(mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
                                 S1.mort.rateratio, time, calibration_data, 
                                 Mortality) + ratediff)
    } 
  }


# Use Mayeda et al. 2017 paper to assign relative mortality of dementia patients
# The mortality rate is HR for mortality among dementa vs. no dementia (whites)
# times mortality rate calculated from lifetables 

mortality.rate.C0D1 <- function(immortal, diff.dem.mort, S1.dem.rateratio, 
                                Sprev, S1.mort.rateratio, time, 
                                calibration_data, Mortality, AD_HR_obs){
  if (immortal == 1){return(0)}
  if (diff.dem.mort == 0){
    return(mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
                               S1.mort.rateratio, time, calibration_data, 
                               Mortality))
  } else{

    #Effect of AD on mortality
    ADHR_t<-Re(polyval(AD_HR_obs$p1,time) / polyval(AD_HR_obs$p2,time))

    return(mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
                               S1.mort.rateratio, time, calibration_data, 
                               Mortality)*ADHR_t)
  }
}

#Combine SEER and Mayeda data to assign mortality rate for 
#cancer AND dementia patients, assuming multiplicative relationship.
mortality.rate.C1D1 <- function(immortal, diff.cancer.mort, diff.dem.mort, 
                                S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                                time, calibration_data, Mortality, AD_HR_obs, 
                                Rel_surv, C1D1mortmultiplier){
  if (immortal == 1){return(0)}
  if (diff.cancer.mort == 0 & diff.dem.mort == 0){
    return(mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
                               S1.mort.rateratio, time, calibration_data, 
                               Mortality))
  }else if (diff.cancer.mort == 0){
    return(mortality.rate.C0D1(immortal, diff.dem.mort, S1.dem.rateratio, Sprev, 
                               S1.mort.rateratio, time, calibration_data, 
                               Mortality, AD_HR_obs))
  }else if (diff.dem.mort == 0){
    return(mortality.rate.C1D0(immortal, diff.cancer.mort, S1.dem.rateratio, 
                               Sprev, S1.mort.rateratio, time, calibration_data, 
                               Mortality, Rel_surv))
  } else{
    
    #Effect of AD on mortality
    ADHR_t<-Re(polyval(AD_HR_obs$p1,time) / polyval(AD_HR_obs$p2,time))
    
    #Effect of cancer on mortality
    rel.surv.5yrs<- Re(polyval(Rel_surv$p1,time) / polyval(Rel_surv$p2,time))
    ratediff<- as.numeric(- 1/5*log(rel.surv.5yrs))
    ratediff<-ifelse(ratediff<0,0,ratediff)
    
    #Multiplier<-Re(polyval(C1D1mortmultiplier$p1,time) / polyval(C1D1mortmultiplier$p2,time))
    
      return(max(((mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev,
                                 S1.mort.rateratio, time, calibration_data,
                                 Mortality)+ratediff)*
                        ADHR_t*C1D1mortmultiplier),
                 mortality.rate.C0D1(immortal, diff.dem.mort, S1.dem.rateratio, Sprev,
                        S1.mort.rateratio, time, calibration_data,
                        Mortality, AD_HR_obs),
                 mortality.rate.C1D0(immortal, diff.cancer.mort, S1.dem.rateratio,
                        Sprev, S1.mort.rateratio, time, calibration_data,
                        Mortality, Rel_surv)))
    
    # return((mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
    #                            S1.mort.rateratio, time, calibration_data, 
    #                            Mortality)+ratediff)*
    #                            ADHR_t*C1D1mortmultiplier)#Multiplier)
    } 
  }


#------------------------------------------------------------------
# Defining mortality rates #
# There are 8 mortality rates for the 8 states 
# States are combos of Cancer, Dementia, and Survival characteristic
#------------------------------------------------------------------

#Define final mortality rates for all 8 states          

# For cancer-free individuals, mortality rate is rate for their dementia 
# status, regardless of S
mortality.rate.C0D0S0 <- function(immortal, S1.dem.rateratio, Sprev, 
                                  S1.mort.rateratio, time, calibration_data, 
                                  Mortality){
  return(mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
                             S1.mort.rateratio, time, calibration_data, 
                             Mortality))}

mortality.rate.C0D0S1 <- function(immortal, S1.dem.rateratio, Sprev, 
                                  S1.mort.rateratio, time, calibration_data, 
                                  Mortality){
  return(mortality.rate.C0D0(immortal, S1.dem.rateratio, Sprev, 
                             S1.mort.rateratio, time, calibration_data, 
                             Mortality))}

mortality.rate.C0D1S0 <- function(immortal, diff.dem.mort, S1.dem.rateratio, 
                                  Sprev, S1.mort.rateratio, time, 
                                  calibration_data, Mortality, AD_HR_obs){
  return(mortality.rate.C0D1(immortal, diff.dem.mort, S1.dem.rateratio, Sprev, 
                             S1.mort.rateratio, time, calibration_data, 
                             Mortality, AD_HR_obs))}

mortality.rate.C0D1S1 <- function(immortal, diff.dem.mort, S1.dem.rateratio, 
                                  Sprev, S1.mort.rateratio, time, 
                                  calibration_data, Mortality, AD_HR_obs){
  return(mortality.rate.C0D1(immortal, diff.dem.mort, S1.dem.rateratio, Sprev, 
                             S1.mort.rateratio, time, calibration_data, 
                             Mortality, AD_HR_obs))}

# For those with cancer and no S, mortality rate is rate for their 
# cancer/dementia status 
mortality.rate.C1D1S0 <- function(immortal, diff.cancer.mort, diff.dem.mort, 
                                  S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                                  time, calibration_data, Mortality, AD_HR_obs, 
                                  Rel_surv,C1D1mortmultiplier){
  return(mortality.rate.C1D1(immortal, diff.cancer.mort, diff.dem.mort, 
                             S1.dem.rateratio, Sprev, S1.mort.rateratio, time, 
                             calibration_data, Mortality, AD_HR_obs, Rel_surv, C1D1mortmultiplier))}

mortality.rate.C1D0S0 <- function(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                  Sprev, S1.mort.rateratio, time, 
                                  calibration_data, Mortality, Rel_surv){
  return(mortality.rate.C1D0(immortal, diff.cancer.mort, S1.dem.rateratio, 
                             Sprev, S1.mort.rateratio, time, calibration_data, 
                             Mortality, Rel_surv))}

# For those with cancer and S = 1, mortality rate is rate for their 
# cancer/dementia status*protective effect of S
mortality.rate.C1D1S1 <- function(immortal, diff.cancer.mort, diff.dem.mort, 
                                  S1.dem.rateratio, Sprev, time, 
                                  calibration_data, Mortality, AD_HR_obs, 
                                  Rel_surv, S1.mort.rateratio, C1D1mortmultiplier){
  return(mortality.rate.C1D1(immortal, diff.cancer.mort, diff.dem.mort, 
                      S1.dem.rateratio, Sprev, S1.mort.rateratio, time, 
                      calibration_data, Mortality, AD_HR_obs, Rel_surv, C1D1mortmultiplier)*
           S1.mort.rateratio)}

mortality.rate.C1D0S1 <- function(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                  Sprev, time, calibration_data, Mortality, 
                                  Rel_surv, S1.mort.rateratio){
  return(mortality.rate.C1D0(immortal, diff.cancer.mort, S1.dem.rateratio, 
                             Sprev, S1.mort.rateratio, time, calibration_data, 
                             Mortality, Rel_surv)*S1.mort.rateratio)}




