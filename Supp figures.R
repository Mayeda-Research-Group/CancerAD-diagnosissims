#Supplementary file figures showing polynomials and calibration:
  


#---- Package Loading ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "ggplot2", "magrittr", "foreign", "deSolve", 
       "numDeriv", "pracma", "knitr")

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


modresults<-compartmodel(immortal = 0,
                         diff.cancer.mort = 1,
                         diff.dem.mort = 1,
                         S1.dem.rateratio = 1,
                         S1.mort.rateratio = 1,
                         C1D1mortmultiplier = .5,
                         Incr_contact=1,
                         Decr_demdx=1,
                         Demdxwaittime=1.5,
                         Sprev = 0,
                         times = seq(from = 0, to = 30, by = 1/12),
                         type = "all")

#Format results
out.data <- modresults$out.data
rate.results <- data.frame(modresults$rate.results)
type<-"all"

if(type == "all" | type == "lung"){
  sex = "all"
} else if(type == "breast"){
  sex = "F"
} else {
  sex = "M"
}

#bring in polynomial calibrations
    
    #Cancer incidence rates -- per Rebecca Graff citations, ok for it to decline over time
    SEER_in_formod<-SEER_inc_obs[,c("time",paste0("Rate_",type))]
    SEER_in_formod$Rate<-SEER_in_formod[,2]
    SEER_in_formod$time<-(SEER_in_formod$time-2.5)
    SEERpoly<-rationalfit(SEER_in_formod$time, SEER_in_formod$Rate, d1=2,d2=0)
    
    #Dementia incidence rates
    AD_in_formod<-AD_inc_obs[,c("Age",paste0("Rate_",sex))]
    AD_in_formod$Rate<-AD_in_formod[,2]
    AD_in_formod$Age<-(AD_in_formod$Age-67.5)
    ADpoly<-rationalfit(AD_in_formod$Age[2:9], AD_in_formod$Rate[2:9], d1=3,d2=1)

    #Cancer survival
    SEER_rel <- SEER_rel.surv[, c("Age", paste0("Rel.surv_", type))] 
    SEER_relsurv_mod<-data.frame(Age=seq(2.5,37.5,5),
                                 Rel.surv=c(rep(SEER_rel[SEER_rel$Age=="65-74",2],2),
                                            rep(SEER_rel[SEER_rel$Age=="75+",2],6)))
    Rel_surv<-rationalfit(SEER_relsurv_mod$Age, SEER_relsurv_mod$Rel.surv, d1=3,d2=0)
    rel.surv.5yrs<- Re(polyval(Rel_surv$p1,seq(0,30,1/12)) / polyval(Rel_surv$p2,seq(0,30,1/12)))
    ratediff<- as.numeric(- 1/5*log(rel.surv.5yrs))
    ratediff<-ifelse(ratediff<0,0,ratediff)
    
    #Dementia survival
    ADHR_formod<-AD_HR_obs[,c("Age","Dementia.HR")]
    ADHR_formod$Age<-c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5)
    AD_HR_obs2<-rationalfit(ADHR_formod$Age, ADHR_formod$Dementia.HR, d1=1,d2=1)

    
#Plot calibrations
  #Cumulative survival
    survplot<-ggplot(out.data, aes(x=time+65,y=(1-DEAD)*100, color="black"))+
      geom_line(size=2)+
      geom_point(data=Lifetables, aes(x=as.numeric(Age), 
                                      y=as.numeric(Cum_surv_obs_all)*100,
                                      color="red"), shape=15, size=4)+
      scale_x_continuous(limits=c(65,95), breaks=seq(65,95,5), minor_breaks = NULL)+
      scale_color_manual(name=NULL, values=c("black", "red"), labels=c("Simulated data", "Calibration data"))+
      ylab("Percent of cohort surviving")+xlab("Age (years)")+
      theme_bw()+
      theme(axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=16), 
            legend.position = c(0.8,0.56),
            legend.text = element_text(size=16))
    survplot
  
  #Dementia incidence rates
    newdem<-out.data$DEMENTIA-lag(out.data$DEMENTIA)
    PTatriskdem<-out.data$PT_C0D0S0+out.data$PT_C0D0S1+out.data$PT_C1D0S0+out.data$C1D0S1
    newPTatriskdem<-PTatriskdem-lag(PTatriskdem)
    out.data$ratedem<-newdem/newPTatriskdem
    
    dem_calibdata<-data.frame(time=seq(0,30,5), demrate_cal=polyval(ADpoly$p1,seq(0,30,5)) / polyval(ADpoly$p2,seq(0,30,5)))
    
    demplot<-ggplot(out.data, aes(x=time+65, y=ratedem*1000, color="black"))+
      geom_line(size=2)+
      geom_point(data=dem_calibdata, aes(x=time+65,y=demrate_cal*1000, color="red"), shape=15, size=4)+
      scale_x_continuous(limits=c(65,95), breaks=seq(65,95,5), minor_breaks = NULL)+
      scale_color_manual(name=NULL, values=c("black", "red"), labels=c("Simulated data", "Calibration data"))+
      ylab("Dementia incidence rate (cases per 1,000 person-years)")+xlab("Age (years)")+
      theme_bw()+
      theme(axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=16), 
            legend.position = c(0.85,0.3),
            legend.text = element_text(size=16))
    demplot
  
  #Cancer incidence rates
    newcan<-out.data$CANCER-lag(out.data$CANCER)
    PTatriskcan<-out.data$PT_C0D0S0+out.data$PT_C0D0S1+out.data$PT_C0D1S0+out.data$C0D1S1+out.data$PT_C0D2S0+out.data$C0D2S1
    newPTatriskcan<-PTatriskcan-lag(PTatriskcan)
    out.data$ratecan<-newcan/newPTatriskcan
    
    can_calibdata<-data.frame(time=seq(0,30,5), canrate_cal=polyval(SEERpoly$p1,seq(0,30,5)) / polyval(SEERpoly$p2,seq(0,30,5)))
    
    canplot<-ggplot(out.data, aes(x=time+65, y=ratecan*1000, color="black"))+
      geom_line(size=2)+
      geom_point(data=can_calibdata, aes(x=time+65,y=canrate_cal*1000, color="red"), shape=15, size=4)+
      scale_x_continuous(limits=c(65,95), breaks=seq(65,95,5), minor_breaks = NULL)+
      scale_color_manual(name=NULL, values=c("black", "red"), labels=c("Simulated data", "Calibration data"))+
      ylab("Cancer incidence rate\n (cases per 1,000 person-years, all cancers combined)")+xlab("Age (years)")+
      theme_bw()+
      theme(axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=16), 
            legend.position = c(0.8,0.3),
            legend.text = element_text(size=16))
  
    canplot
    
    #Dementia mortality HR
    rate.nodem.mort <- #weighted average rate of death from no dementia states
      (rate.results$dDEADC0D0S0)/(out.data$C0D0S0)
    rate.dem.mort <- #weighted average rate of death from no dementia states
      (rate.results$dDEADC0D1S0)/(out.data$C0D1S0)
    rate.results$irrdemort_sim<-rate.dem.mort/rate.nodem.mort
    
    adhr_calibdata<-data.frame(time=seq(0,30,5), adhr_cal=polyval(AD_HR_obs2$p1,seq(0,30,5)) / polyval(AD_HR_obs2$p2,seq(0,30,5)))
    
    demHRplot<-ggplot(rate.results, aes(x=time+65, y=irrdemort_sim, color="black"))+
      geom_line(size=2)+
      geom_point(data=adhr_calibdata, aes(x=time+65,y=adhr_cal, color="red"), shape=15, size=4)+
      scale_x_continuous(limits=c(65,95), breaks=seq(65,95,5), minor_breaks = NULL)+
      scale_color_manual(name=NULL, values=c("black", "red"), labels=c("Simulated data", "Calibration data"))+
      ylab("Hazard ratio for dementia on mortality")+xlab("Age (years)")+
      theme_bw()+
      theme(axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=16), 
            legend.position = c(0.8,0.5),
            legend.text = element_text(size=16))
    
    demHRplot  
  
#Check cancer mortality rate difference
    
    rate.nocan.mort <- #weighted average rate of death from no cancer
      (rate.results$dDEADC0D0S0)/(out.data$C0D0S0)
    rate.can.mort <- #weighted average rate of death from cancer 
      (rate.results$dDEADC1D0S0)/(out.data$C1D0S0)
    rate.results$ratediff_sim<-rate.can.mort-rate.nocan.mort
    
    canRD_calibdata<-data.frame(time=seq(0,30,5), canRD_cal=ratediff[rate.results$time %in% seq(0,30,5)])
    
    canRDplot<-ggplot(rate.results, aes(x=time+65, y=ratediff_sim, color="black"))+
      geom_line(size=2)+
      geom_point(data=canRD_calibdata, aes(x=time+65,y=canRD_cal, color="red"), shape=15, size=4)+
      scale_x_continuous(limits=c(65,95), breaks=seq(65,95,5), minor_breaks = NULL)+
      scale_color_manual(name=NULL, values=c("black", "red"), labels=c("Simulated data", "Calibration data"))+
      ylab("Rate difference for cancer on mortality")+xlab("Age (years)")+
      theme_bw()+
      theme(axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=16), 
            legend.position = c(0.8,0.5),
            legend.text = element_text(size=16))
    
    canRDplot
    
    ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Supp/survplot.jpg", 
           device="jpg", plot=survplot,  dpi="retina")
    
    ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Supp/demplot.jpg", 
           device="jpg", plot=demplot,  dpi="retina")
    
    ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Supp/canplot.jpg", 
           device="jpg", plot=canplot,  dpi="retina")
    
    ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Supp/demHRplot.jpg", 
           device="jpg", plot=demHRplot,  dpi="retina")
    
    ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Supp/canRDplot.jpg", 
           device="jpg", plot=canRDplot,  dpi="retina")
    