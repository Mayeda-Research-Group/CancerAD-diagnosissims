#------------------------------------------------------------------
# Title: Run simulations
# Author: Eleanor Hayes-Larson
# Co-author/editor: Crystal Shaw
#------------------------------------------------------------------

#---- Package Loading ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "ggplot2", "magrittr", "foreign", "deSolve", 
       "numDeriv", "pracma")

#------------------------------------------------------------------
#---- Load model ----
#------------------------------------------------------------------
source(here("Code", "Multistate model", "Rate inputs.R"))
#Script that generates initial states
source(here("Code", "Multistate model", "gen_init.R")) 
#Script for cancer incidence rates
source(here("Code", "Multistate model", "get_cancer.R"))
#Script for dementia incidence rates
source(here("Code", "Multistate model", "get_dem.R"))
#Script for dementia diagnosis rates
source(here("Code", "Multistate model", "get_demdx.R"))
#Script for mortality rates
source(here("Code", "Multistate model", "mortality_rates.R"))
#Script for dementia incidence rates
source(here("Code", "Multistate model", "diff_eq_model.R"))
source(here("Code", "Multistate model", "Model.R"))

#------------------------------------------------------------------
# Create function that runs model and produces results. 
#------------------------------------------------------------------
multirun <- function( 
  arg_S1.dem.rateratio, arg_S1.mort.rateratio, arg_Sprev, arg_C1D1mortmultiplier,
  arg_Incr_contact, 
  arg_Decr_demdx, arg_Demdxwaittime,
  cancertype){
  modresults <- list()
  
  #Run model
  modresults <- 
    compartmodel(immortal = 0, 
                 diff.cancer.mort = 1, 
                 diff.dem.mort = 1, 
                 S1.dem.rateratio = arg_S1.dem.rateratio, 
                 S1.mort.rateratio = arg_S1.mort.rateratio, 
                 C1D1mortmultiplier = arg_C1D1mortmultiplier,
                 Incr_contact=arg_Incr_contact, 
                 Decr_demdx=arg_Decr_demdx,
                 Demdxwaittime = arg_Demdxwaittime,
                 Sprev = arg_Sprev,
                 times = seq(from = 0, to = 30, by = 1/12),
                 type = cancertype)
  
  #Format results
  out.data <- modresults$out.data
  rate.results <- data.frame(modresults$rate.results)
  
  #Get model results for any dementia (undiagnosed)
  #Get lifetime relative risk
  relativerisk <- 
    data.frame(time = out.data$time, 
               RR = ((out.data$DEMWC/out.data$CBEFOREDEM)/
                       (out.data$DEMNOC/(1-out.data$CBEFOREDEM))))
  lifetimeRR <- relativerisk$RR[relativerisk$time == 30] #request RR for very large time to ensure steady state--i.e. all have died.
  lifetimeRR
  
  #Get time-specific rate ratio from diffEQ model output
  IRR.months <- rate.results$dDEMWC/
    rate.results$dDEMNOC*
    ((out.data$C0D0S0+out.data$C0D0S1)/
       (out.data$C1D0S0+out.data$C1D0S1))
  weights <- out.data$C0D0S0 + out.data$C1D0S0 + out.data$C0D0S1 + 
    out.data$C1D0S1
  
  overallIRR <- exp(weighted.mean(log(IRR.months[2:361]), weights[2:361]))
  overallIRR
  
  decade.IRRs <- rep(NA, 3)
  for (i in 1:3){
    if (i == 1){
      decade.IRRs[i] <- exp(weighted.mean(log(IRR.months[2:(i*120)]), 
                                          weights[2:(i*120)]))}
    else{
      decade.IRRs[i] <- 
        exp(weighted.mean(log(IRR.months[((i-1)*120 + 1):(i*120)]), 
                          weights[((i-1)*120 + 1):(i*120)]))
    }
  }
  decade.IRRs
  
  #rate calculated from PT
  #Overall
  out.data$wC_n_mo<-out.data$DEMWC-lag(out.data$DEMWC)
  
  out.data$noC_n_mo<-out.data$DEMNOC-lag(out.data$DEMNOC)
  
  
  out.data$wC_PT_mo<-(out.data$PT_C1D0S0+out.data$PT_C1D0S1)-
    (lag(out.data$PT_C1D0S0)+lag(out.data$PT_C1D0S1))
  
  out.data$noC_PT_mo<-(out.data$PT_C0D0S0+out.data$PT_C0D0S1)-
    (lag(out.data$PT_C0D0S0)+lag(out.data$PT_C0D0S1))
  
  
  out.data$Monthly_IRR_PTcalc<-(out.data$wC_n_mo/out.data$wC_PT_mo)/
    (out.data$noC_n_mo/out.data$noC_PT_mo)
  
  Monthly_IRR_PTcalc<-out.data$Monthly_IRR_PTcalc
  
  IRR_PTcalc <- exp(weighted.mean(log(Monthly_IRR_PTcalc[2:361]), weights[2:361]))
  IRR_PTcalc
  
  Decade_IRR_PTcalc <- rep(NA, 3)
  
  for (i in 1:3){
    
    if (i == 1){
      Decade_IRR_PTcalc[i] <- exp(weighted.mean(log(Monthly_IRR_PTcalc[2:(i*120)]), 
                                                weights[2:(i*120)]))}
    else{
      Decade_IRR_PTcalc[i] <- 
        exp(weighted.mean(log(Monthly_IRR_PTcalc[((i-1)*120 + 1):(i*120)]), 
                          weights[((i-1)*120 + 1):(i*120)]))
    }
  }
  Decade_IRR_PTcalc
  
  
  #Obtain results for diagnosed dementia:
  
  #Get time-specific rate ratio from diffEQ model output conditional on undxed dem
  Demdx_IRR.months <- rate.results$dDEMdxWC/
    rate.results$dDEMdxNOC*
    ((out.data$C0D1S0+out.data$C0D1S1)/
       (out.data$C1D1S0+out.data$C1D1S1))
  Demdx_weights <- out.data$C0D0S0 + out.data$C1D0S0 + out.data$C0D0S1 + 
    out.data$C1D0S1 + out.data$C0D1S0 + out.data$C1D1S0 + out.data$C0D1S1 + 
    out.data$C1D1S1
  
  Demdx_overallIRR <- exp(weighted.mean(log(Demdx_IRR.months[2:361]), Demdx_weights[2:361]))
  Demdx_overallIRR
  
  Demdx_decade.IRRs <- rep(NA, 3)
  for (i in 1:3){
    if (i == 1){
      Demdx_decade.IRRs[i] <- exp(weighted.mean(log(Demdx_IRR.months[2:(i*120)]), 
                                                Demdx_weights[2:(i*120)]))}
    else{
      Demdx_decade.IRRs[i] <- 
        exp(weighted.mean(log(Demdx_IRR.months[((i-1)*120 + 1):(i*120)]), 
                          Demdx_weights[((i-1)*120 + 1):(i*120)]))
    }
  }
  Demdx_decade.IRRs
  
  
  
  #PT calcs
  
  out.data$Demdx_wc_n_mo<-out.data$DEMdxWC-lag(out.data$DEMdxWC)
  out.data$Demdx_noC_n_mo<-out.data$DEMdxNOC-lag(out.data$DEMdxNOC)
  
  
  #Conditional on undxed dem
  out.data$Demdx_wc_PTwithD_mo<-(out.data$PT_C1D1S0+out.data$PT_C1D1S1)-
    (lag(out.data$PT_C1D1S0)+lag(out.data$PT_C1D1S1))
  out.data$Demdx_noC_PTwithD_mo<-(out.data$PT_C0D1S0+out.data$PT_C0D1S1)-
    (lag(out.data$PT_C0D1S0)+lag(out.data$PT_C0D1S1))
  
  out.data$Demdx_Monthly_IRR_wD_PTcalc<-(out.data$Demdx_wc_n_mo/out.data$Demdx_wc_PTwithD_mo)/
    (out.data$Demdx_noC_n_mo/out.data$Demdx_noC_PTwithD_mo)
  
  Demdx_Monthly_IRR_wD_PTcalc<-out.data$Demdx_Monthly_IRR_wD_PTcalc
  
  Demdx_Monthly_IRR_wD_PTcalc<-exp(weighted.mean(log(Demdx_Monthly_IRR_wD_PTcalc[2:361]), Demdx_weights[2:361]))
  Demdx_Monthly_IRR_wD_PTcalc
  
  #Unconditional on dementia
  
  out.data$Demdx_wc_PT_mo<-(out.data$PT_C1D0S0+out.data$PT_C1D0S1+out.data$PT_C1D1S0+out.data$PT_C1D1S1)-
    (lag(out.data$PT_C1D0S0)+lag(out.data$PT_C1D0S1)+lag(out.data$PT_C1D1S0)+lag(out.data$PT_C1D1S1))
  out.data$Demdx_noC_PT_mo<-(out.data$PT_C0D0S0+out.data$PT_C0D0S1+out.data$PT_C0D1S0+out.data$PT_C0D1S1)-
    (lag(out.data$PT_C0D0S0)+lag(out.data$PT_C0D0S1)+lag(out.data$PT_C0D1S0)+lag(out.data$PT_C0D1S1))
  
  out.data$Demdx_Monthly_IRR_PTcalc<-(out.data$Demdx_wc_n_mo/out.data$Demdx_wc_PT_mo)/
    (out.data$Demdx_noC_n_mo/out.data$Demdx_noC_PT_mo)
  
  Demdx_Monthly_IRR_PTcalc<-out.data$Demdx_Monthly_IRR_PTcalc
  #plot(Demdx_Monthly_IRR_PTcalc)
  Demdx_IRR_PTcalc <- exp(weighted.mean(log(Demdx_Monthly_IRR_PTcalc[2:361]), Demdx_weights[2:361]))
  Demdx_IRR_PTcalc
  
  Demdx_decade_IRR_PTcalc <- rep(NA, 3)
  
  for (i in 1:3){
    
    if (i == 1){
      Demdx_decade_IRR_PTcalc[i] <- exp(weighted.mean(log(Demdx_Monthly_IRR_PTcalc[2:(i*120)]), 
                                                      Demdx_weights[2:(i*120)]))}
    else{
      Demdx_decade_IRR_PTcalc[i] <- 
        exp(weighted.mean(log(Demdx_Monthly_IRR_PTcalc[((i-1)*120 + 1):(i*120)]), 
                          Demdx_weights[((i-1)*120 + 1):(i*120)]))
    }
  }
  Demdx_decade_IRR_PTcalc
  
  
  
  #Return model parameter values and results
  return(list(cancertype = paste(cancertype), 
              S1.dem.rateratio = arg_S1.dem.rateratio, 
              S1.mort.rateratio = arg_S1.mort.rateratio, 
              Sprev = arg_Sprev, 
              Incr_contact = arg_Incr_contact, 
              Decr_demdx = arg_Decr_demdx, 
              Demdxwaittime = arg_Demdxwaittime,
              Mortmultfactor = arg_C1D1mortmultiplier,
              
              lifetimeRR = lifetimeRR,
              IRR.overall = overallIRR,
              IRR.6575 = decade.IRRs[1],
              IRR.7585 = decade.IRRs[2], 
              IRR.8595 = decade.IRRs[3], 

              IRR_PTcalc = IRR_PTcalc,
              IRR.6575_PTcalc = Decade_IRR_PTcalc[1],
              IRR.7585_PTcalc = Decade_IRR_PTcalc[2], 
              IRR.8595_PTcalc = Decade_IRR_PTcalc[3], 

              Demdx_overallIRR = Demdx_overallIRR,
              Demdx_IRR.6575 = Demdx_decade.IRRs[1],
              Demdx_IRR.7585 = Demdx_decade.IRRs[2], 
              Demdx_IRR.8595 = Demdx_decade.IRRs[3], 

              Demdx_Monthly_IRR_wD_PTcalc=Demdx_Monthly_IRR_wD_PTcalc,
              
              Demdx_IRR_PTcalc=Demdx_IRR_PTcalc,
              Demdx_IRR.6575_PTcalc = Demdx_decade_IRR_PTcalc[1],
              Demdx_IRR.7585_PTcalc = Demdx_decade_IRR_PTcalc[2], 
              Demdx_IRR.8595_PTcalc = Demdx_decade_IRR_PTcalc[3], 

                            outdata = out.data,
              rate.results = rate.results))
}

#------------------------------------------------------------------
# Define parameter values of interest 
#------------------------------------------------------------------
S1demvals <- c(1)#,0.7)
S1mortvals <- c(1)#,0.7)
Sprevvals <- c(0)#,0.3)
Mortmult_vals<-c(0.7, 0.6, 0.5, 0.4)
typevals <- c("all","lung", "breast", "prostate")
Incr_contactvals<-c(1, 1.1, 1.2, 1.3, 1.4, 1.5)
Decr_demdxvals<-c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
Demdxwaittimevals<-c(.05, 0.5, 1, 1.5, 2)
args <- expand.grid(S1demvals = S1demvals, S1mortvals = S1mortvals, 
                    Sprevvals = Sprevvals, Mortmult_vals = Mortmult_vals, 
                    typevals = typevals, 
                    Incr_contactvals=Incr_contactvals, Decr_demdxvals=Decr_demdxvals,
                    Demdxwaittimevals=Demdxwaittimevals)
args<-args[((args$S1demvals==1 & args$S1mortvals==1 & args$Sprevvals==0 )|
             (args$S1demvals==0.7 & args$S1mortvals==0.7 & args$Sprevvals==0.3 )) &
             (args$Incr_contactvals==1 | args$Decr_demdxvals==1 ),]
#------------------------------------------------------------------
# Run model across parameter values and store results 
#------------------------------------------------------------------
#Timing the code
#start <- Sys.time()
Allsims_dxbias <- mapply(FUN = multirun, 
                         arg_S1.dem.rateratio = args$S1demvals, 
                         arg_S1.mort.rateratio = args$S1mortvals, 
                         arg_Sprev = args$Sprevvals, 
                         arg_C1D1mortmultiplier = args$Mortmult_vals,
                         arg_Incr_contact = args$Incr_contactvals,
                         arg_Decr_demdx = args$Decr_demdxvals,
                         arg_Demdxwaittime=args$Demdxwaittimevals,
                         cancertype = args$typevals)




#End code timing
#end <- Sys.time() 
#end-start

save(Allsims_dxbias, file = here("Output", "Allsims_dxbias.Rdata"))

