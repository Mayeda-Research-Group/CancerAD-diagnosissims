#------------------------------------------------------------------
# Part 3: Specify model
#      #Define the differential equations for the model#
#      #Define initial conditions for the model#
#------------------------------------------------------------------
# diff.eq.model() is a function which defines the set of differential equations
# inputs: time, states, parameters
## time: a vector of times at which you want the derivates
## states: a vector of number of people in each state with the entries labeled
## parameters: a list of parameters that other functions in diff.eq.model needs
### Cancer_inc, time, S1.dem.rateratio, Sprev, calibration_data, immortal, 
### Mortality, diff.cancer.mort, diff.dem.mort, AD_HR_obs, Rel_surv, 
### S1.mort.rateratio

# outputs a list where the first entry is a vector derivatives and the second entry is NULL
## this output format is required by the differential equation solver

diff.eq.model<- function(time, states, parms){
  with(as.list(c(states, parms)),
       {
         #Rates for 9 states in the model
         dC0D0S0 <- -get.cancer(time, Cancer_inc)*C0D0S0 - 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C0D0S0 - 
           mortality.rate.C0D0S0(immortal, S1.dem.rateratio, Sprev, 
                                 S1.mort.rateratio, time, calibration_data, 
                                 Mortality)*C0D0S0
         
         dC0D1S0 <- -get.cancer(time, Cancer_inc)*C0D1S0 + 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C0D0S0 - 
           get.demdx.C0(time, Demdx_inc)*C0D1S0 -
           mortality.rate.C0D1S0(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, AD_HR_obs)*C0D1S0
         
         dC1D0S0 <-  get.cancer(time, Cancer_inc)*C0D0S0 - 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C1D0S0 - 
           mortality.rate.C1D0S0(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, Rel_surv)*C1D0S0
         
         dC1D1S0 <-  get.cancer(time, Cancer_inc)*C0D1S0 + 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C1D0S0 - 
           get.demdx.C1(time, Incr_contact, Decr_demdx, Demdx_inc)*C1D1S0 -
           mortality.rate.C1D1S0(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                                 time, calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, C1D1mortmultiplier)*C1D1S0
         
         dC0D0S1 <- -get.cancer(time, Cancer_inc)*C0D0S1 - 
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C0D0S1 - 
           mortality.rate.C0D0S1(immortal, S1.dem.rateratio, Sprev, 
                                 S1.mort.rateratio, time, calibration_data, 
                                 Mortality)*C0D0S1
         
         dC0D1S1 <- -get.cancer(time, Cancer_inc)*C0D1S1 + 
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C0D0S1 - 
           get.demdx.C0(time, Demdx_inc)*C0D1S1 -
           mortality.rate.C0D1S1(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, AD_HR_obs)*C0D1S1
         
         dC1D0S1 <-  get.cancer(time, Cancer_inc)*C0D0S1 - 
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C1D0S1 - 
           mortality.rate.C1D0S1(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                 Sprev, time, calibration_data, Mortality, 
                                 Rel_surv, S1.mort.rateratio)*C1D0S1
         
         dC1D1S1 <-  get.cancer(time, Cancer_inc)*C0D1S1 + 
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C1D0S1 -
           get.demdx.C1(time, Incr_contact, Decr_demdx, Demdx_inc)*C1D1S1 -
           mortality.rate.C1D1S1(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, time, 
                                 calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, S1.mort.rateratio, C1D1mortmultiplier)*C1D1S1
         
         #Diagnosed dementia states
         dC0D2S0 <- -get.cancer(time, Cancer_inc)*C0D2S0 + 
           get.demdx.C0(time, Demdx_inc)*C0D1S0 - 
           mortality.rate.C0D1S0(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, AD_HR_obs)*C0D2S0
         
         dC1D2S0 <- get.cancer(time, Cancer_inc)*C0D2S0 +
           get.demdx.C1(time, Incr_contact, Decr_demdx, 
                        Demdx_inc)*C1D1S0 - 
           mortality.rate.C1D1S0(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                                 time, calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, C1D1mortmultiplier)*C1D2S0
         
         dC0D2S1 <- -get.cancer(time, Cancer_inc)*C0D2S1 + 
           get.demdx.C0(time, Demdx_inc)*C0D1S1 - 
           mortality.rate.C0D1S1(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, AD_HR_obs)*C0D2S1
         
         dC1D2S1 <- get.cancer(time, Cancer_inc)*C0D2S1 + 
           get.demdx.C1(time, Incr_contact, Decr_demdx, 
                        Demdx_inc)*C1D1S1 - 
           mortality.rate.C1D1S1(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, time, 
                                 calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, S1.mort.rateratio, C1D1mortmultiplier)*C1D2S1
         
         dDEAD <-  mortality.rate.C0D0S0(immortal, S1.dem.rateratio, Sprev, 
                                         S1.mort.rateratio, time, 
                                         calibration_data, Mortality)*C0D0S0 + 
           mortality.rate.C1D0S0(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, 
                                 Rel_surv)*C1D0S0 + 
           mortality.rate.C0D1S0(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, 
                                 AD_HR_obs)*(C0D1S0 + C0D2S0)+ 
           mortality.rate.C1D1S0(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                                 time, calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, C1D1mortmultiplier)*(C1D1S0 + C1D2S0) +
           mortality.rate.C0D0S1(immortal, S1.dem.rateratio, Sprev, 
                                 S1.mort.rateratio, time, calibration_data, 
                                 Mortality)*C0D0S1 + 
           mortality.rate.C1D0S1(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                 Sprev, time, calibration_data, Mortality, 
                                 Rel_surv, S1.mort.rateratio)*C1D0S1 + 
           mortality.rate.C0D1S1(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, 
                                 AD_HR_obs)*(C0D1S1 + C0D2S1) + 
           mortality.rate.C1D1S1(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, time, 
                                 calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, S1.mort.rateratio, C1D1mortmultiplier)*(C1D1S1 + C1D2S1)
         
         
         #Person time tallies for each state:
         dPT_C0D0S0<-C0D0S0
         dPT_C0D0S1<-C0D0S1
         dPT_C1D0S0<-C1D0S0
         dPT_C1D0S1<-C1D0S1
         dPT_C0D1S0<-C0D1S0
         dPT_C0D1S1<-C0D1S1
         dPT_C1D1S0<-C1D1S0
         dPT_C1D1S1<-C1D1S1
         dPT_C0D2S0<-C0D2S0
         dPT_C0D2S1<-C0D2S1
         dPT_C1D2S0<-C1D2S0
         dPT_C1D2S1<-C1D2S1
         
         #The rest of these are for tallying results (e.g. cumul. incidence)
         dDEMNOC <- 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C0D0S0 + 
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C0D0S1 
         
         dDEMWC <- 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C1D0S0 + 
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*C1D0S1
         
         dDEMdxNOC <-
           get.demdx.C0(time, Demdx_inc)*(C0D1S0+C0D1S1)
         
         dDEMdxWC <-
           get.demdx.C1(time, Incr_contact, Decr_demdx,
                        Demdx_inc)*(C1D1S0+C1D1S1)
         
         dDEMENTIA <- 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*(C0D0S0 + C1D0S0) + 
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*(C0D0S1 + C1D0S1)
         
         dDEMENTIAdx <- get.demdx.C0(time, Demdx_inc)*(C0D1S0+C0D1S1) + 
           get.demdx.C1(time, Incr_contact, Decr_demdx,
                        Demdx_inc)*(C1D1S0+C1D1S1)
         
         dCANCER <- 
           get.cancer(time, Cancer_inc)*(C0D0S0 + C0D1S0 + C0D0S1 + C0D1S1 + C0D2S0 + C0D2S1)
         
         dDEADC0D0S0 <- 
           mortality.rate.C0D0S0(immortal, S1.dem.rateratio, Sprev, 
                                 S1.mort.rateratio, time, calibration_data, 
                                 Mortality)*C0D0S0
         
         dDEADC1D0S0 <- 
           mortality.rate.C1D0S0(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, Rel_surv)*C1D0S0
         
         dDEADC0D1S0 <- 
           mortality.rate.C0D1S0(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, AD_HR_obs)*C0D1S0
         
         dDEADC1D1S0 <- 
           mortality.rate.C1D1S0(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                                 time, calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, C1D1mortmultiplier)*C1D1S0
         
         dDEADC0D0S1 <- 
           mortality.rate.C0D0S1(immortal, S1.dem.rateratio, Sprev, 
                                 S1.mort.rateratio, time, calibration_data, 
                                 Mortality)*C0D0S1
         
         dDEADC1D0S1 <- 
           mortality.rate.C1D0S1(immortal, diff.cancer.mort, S1.dem.rateratio, 
                                 Sprev, time, calibration_data, Mortality, 
                                 Rel_surv, S1.mort.rateratio)*C1D0S1
         
         dDEADC0D1S1 <- 
           mortality.rate.C0D1S1(immortal, diff.dem.mort, S1.dem.rateratio, 
                                 Sprev, S1.mort.rateratio, time, 
                                 calibration_data, Mortality, AD_HR_obs)*C0D1S1
         
         dDEADC1D1S1 <- 
           mortality.rate.C1D1S1(immortal, diff.cancer.mort, diff.dem.mort, 
                                 S1.dem.rateratio, Sprev, time, 
                                 calibration_data, Mortality, AD_HR_obs, 
                                 Rel_surv, S1.mort.rateratio, C1D1mortmultiplier)*C1D1S1
         
         dCBEFOREDEM <- get.cancer(time, Cancer_inc)*(C0D0S0 + C0D0S1)
         
         dDEMNOSNOC <- 
           get.dem.S0(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*(C0D0S0) 
         
         dDEMWITHSNOC <-  
           get.dem.S1(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                      calibration_data, Dem_inc)*(C0D0S1)
         
         list(c(dC0D0S0, dC0D1S0, dC1D0S0, dC1D1S0, dC0D0S1, dC0D1S1, dC1D0S1, 
                dC1D1S1, dC0D2S0, dC0D2S1, dC1D2S0, dC1D2S1, 
                dDEAD, 
                dPT_C0D0S0, dPT_C0D0S1, dPT_C1D0S0, dPT_C1D0S1, 
                dPT_C0D1S0, dPT_C0D1S1, dPT_C1D1S0, dPT_C1D1S1, 
                dPT_C0D2S0, dPT_C0D2S1, dPT_C1D2S0, dPT_C1D2S1,
                dDEMNOC, dDEMWC, dDEMdxNOC, dDEMdxWC, dDEMENTIA, dDEMENTIAdx, 
                dCANCER, 
                dDEADC0D0S0, dDEADC1D0S0, dDEADC0D1S0, dDEADC1D1S0, dDEADC0D0S1, 
                dDEADC1D0S1, dDEADC0D1S1, dDEADC1D1S1, dCBEFOREDEM, dDEMNOSNOC, 
                dDEMWITHSNOC), NULL)
       })}
