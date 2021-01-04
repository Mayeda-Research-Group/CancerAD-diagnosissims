#Examining simulation results

#---- Package Loading ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "ggplot2", "magrittr", "foreign", "deSolve", 
       "numDeriv", "pracma", "dplyr", "RColorBrewer")


# Load sim results
load(here("Output", "Allsims_dxbias.Rdata"))


results<-data.frame(t(data.frame(Allsims_dxbias)))
exampledata<-data.frame(results[results$cancertype=="all" & 
                                  results$Incr_contact==1 & 
                                  results$Decr_demdx==1  & 
                                  results$Mortmultfactor==0.5 &
                                  results$Demdxwaittime==1.5,"outdata"])
colnames(exampledata)<-substr(colnames(exampledata),6,10000)  #fixing col names
exampledata$C0D0<-exampledata$C0D0S0+exampledata$C0D0S1
exampledata$C0D1<-exampledata$C0D1S0+exampledata$C0D1S1
exampledata$C0D2<-exampledata$C0D2S0+exampledata$C0D2S1
exampledata$C1D0<-exampledata$C1D0S0+exampledata$C1D0S1
exampledata$C1D1<-exampledata$C1D1S0+exampledata$C1D1S1
exampledata$C1D2<-exampledata$C1D2S0+exampledata$C1D2S1
exampledata$D0<-exampledata$C0D0+exampledata$C1D0
exampledata$D1<-exampledata$C0D1+exampledata$C1D1
exampledata$D2<-exampledata$C0D2+exampledata$C1D2
test2<-c("time","DEMENTIA", "DEMENTIAdx","DEAD")
to.plot2<-gather(exampledata[,test2], "state","number",-time)

#plot cohort over time.
Cohort_buildfinal2<-ggplot(to.plot2,aes(time+65,number,
                                      group=factor(state, levels=c("DEAD", "DEMENTIA", "DEMENTIAdx")),
                                      color=factor(state, levels=c("DEAD", "DEMENTIA", "DEMENTIAdx"))))+
  geom_line(size=2) + ylab("Cumulative proportion of cohort experiencing \ndementia, dementia diagnosis, and death")+theme_bw()+
  scale_x_continuous(breaks=c(65, 75, 85, 95))+
  scale_color_manual(name=NULL,values=c("darkred","lightblue", "blue"), 
                     labels=c("Death", 
                              "Dementia", 
                              "Diagnosed dementia"))+
  xlab("Age")+theme(legend.position = c(0.8,0.5),
                                          legend.text = element_text(size=16),
                                          axis.text.x = element_text(size=20), 
                                          axis.text.y = element_text(size=20), 
                                          axis.title.x = element_text(size=20, face="bold"), 
                                          axis.title.y = element_text(size=20, face="bold"))

Cohort_buildfinal2

#need to add save of this figure:
ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Figure1.jpg", 
       device="jpg", plot=Cohort_buildfinal2, dpi="retina")


#Calculate cumulative proportion dead, dementia, and dementia dx at 30 years
exampledata[exampledata$time==30, c("DEMENTIA", "DEMENTIAdx", "DEAD")]
exampledata[exampledata$time==30,"DEMENTIAdx"]/exampledata[exampledata$time==30,"DEMENTIA"]


#Main results

results<-results[,1:26]
results$dxbias<-as.numeric(results$Incr_contact)*as.numeric(results$Decr_demdx)

#print main results

basecase_nondiff<-results[results$Incr_contact==1 & 
                            results$Decr_demdx==1 & results$Demdxwaittime==1.5 & 
                            results$Mortmultfactor==.5,]

print(basecase_nondiff[,c("cancertype", "Demdx_IRR_PTcalc")])




#Input Ording et al. data
ording_AD<-data.frame(cancertype=c("All", "Lung", "Breast","Prostate"), 
                   Est = c(0.94, 0.84, 0.95, 0.96), 
                   LCI = c(0.92,0.72, 0.91, 0.89), 
                   UCI = c(0.96, 0.97, 1.00,1.03),
                   Data = rep("Ording et al. (AD only)",4)
)

ording_all<-data.frame(cancertype=c("All", "Lung", "Breast","Prostate"), 
                   Est =  c(0.96, 1.12, 1.00, 0.97), 
                   LCI =  c(0.95, 1.04, 0.97, 0.93), 
                   UCI =  c(0.97, 1.22, 1.03,1.01),
                   Data = rep("Ording et al. (All-cause dementia)",4)
)


#Need to check with Monica if I shoudl use All cancer types , or summary HR for all studies
ospina<-data.frame(cancertype=c("All", "Breast","Prostate"), 
                   Est = c(0.81, 0.93, 0.99), 
                   LCI = c(0.70, 0.87, 0.87), 
                   UCI  = c(0.94, 0.98, 1.13),
                   Data = rep("Ospina-Romero et al.",3)
)




#MSE to determine best-fit models

    sim_res_dxbias<-results %>% filter(Demdxwaittime==1.5, Mortmultfactor==0.5) %>% mutate(cancertype=factor(str_to_sentence(cancertype)))
    obs_data<-rbind(ording_AD, ording_all, ospina)
    
    obs_data$cancertype <- factor(obs_data$cancertype)
    
    sim_obs_merge<-full_join(sim_res_dxbias, obs_data, by = "cancertype")
    
    
    sim_obs_merge$sim_obs_diffsq<-(as.numeric(sim_obs_merge$Demdx_IRR_PTcalc)-as.numeric(sim_obs_merge$Est))^2
    
    bestfit<-sim_obs_merge %>% filter(cancertype!="Lung") %>% group_by(cancertype, dxbias) %>% summarize(MSE=weighted.mean(sim_obs_diffsq, c(0.5, 0.5, 1)))
    bestfit2<-sim_obs_merge %>% filter(cancertype=="Lung") %>% group_by(cancertype, dxbias) %>% summarize(MSE=weighted.mean(sim_obs_diffsq, c(0.5, 0.5)))
    bestfit3<-rbind(bestfit, bestfit2)
    
    bestfit4<-bestfit3 %>% group_by(cancertype) %>% filter(MSE==min(MSE))
    
#Format results for plotting
    Nondif_delay<-results %>% filter(Decr_demdx==1 & Incr_contact==1 & 
                                       Mortmultfactor==0.5 & Demdxwaittime == 1.5) %>%
                    mutate(cancertype=str_to_sentence(cancertype), Est=Demdx_IRR_PTcalc, LCI=NA, UCI=NA, Data = "Non-differential Delay simulation") %>%
                    select(cancertype, Est, LCI, UCI, Data)
    
    
    
    Dif_delay<-left_join(bestfit4, sim_res_dxbias, by=c("cancertype", "dxbias"))
    Dif_delay<-Dif_delay %>% mutate(cancertype=str_to_sentence(cancertype),
                                    Est=Demdx_IRR_PTcalc, LCI=NA, UCI=NA, 
                                    Data = "Best-match Differential Delay simulation") %>%
                                    select(cancertype, Est, LCI, UCI, Data)
      
    
    Overlaydata<-rbind(obs_data, Nondif_delay, 
                       Dif_delay)
    Overlaydata$Data=factor(Overlaydata$Data, levels=c("Non-differential Delay simulation",
                                                       "Best-match Differential Delay simulation",
                                                       "Ospina-Romero et al.",
                                                       "Ording et al. (AD only)",
                                                       "Ording et al. (All-cause dementia)"))

#Plot    
    overlay<-ggplot(data=Overlaydata, aes(x=1, 
                                            group=Data, color=Data, shape=Data))+
      geom_pointrange(aes(y=as.numeric(Est), ymin=LCI, ymax=UCI), size=1.5, position=position_dodge(.5))+
      scale_shape_manual(name="",values=c(17, 15, 16, 16, 16))+
      scale_color_manual(name="", values=c("dodgerblue3", "dodgerblue4", "orange", "tomato3", "tomato4"))+
      scale_x_discrete(labels=c("All", "Breast", "Prostate", "Lung"))+
      scale_y_continuous(breaks=c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2))+
      geom_hline(yintercept=1, colour="black", lwd=1) +
      theme_bw()+facet_grid(~factor(unlist(cancertype), 
                                    levels=c("All", "Breast","Prostate", "Lung")))+
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            legend.text = element_text(size=16),
            axis.text.x = element_text(size=20), 
            axis.text.y = element_text(size=20), 
            axis.title.x = element_text(size=20, face="bold"), 
            axis.title.y = element_text(size=20, face="bold"),
            strip.text.x = element_text(size=16))+
      ylab("Effect estimate")+xlab("Cancer type")+
      guides(color=guide_legend(nrow=5), shape=guide_legend(nrow=5))
    
    overlay
    
    ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Figure2.jpg", 
           device="jpg", plot=overlay, dpi="retina")




#Plot sensitivity analyses
nodxbias<-ggplot(data=results[results$cancertype=="all" & results$Incr_contact==1 &results$Decr_demdx==1,])+
  geom_point(aes(x=unlist(Demdxwaittime), y=unlist(Demdx_IRR_PTcalc), 
                 group=factor(unlist(Mortmultfactor)), 
                 color=factor(unlist(Mortmultfactor))),
             size=2)+
  geom_line(aes(x=unlist(Demdxwaittime), y=unlist(Demdx_IRR_PTcalc), 
                group=factor(unlist(Mortmultfactor)), 
                color=factor(unlist(Mortmultfactor))), 
                size=2)+
  ylim(0.7,1.1)+
  scale_x_continuous(breaks=seq(0,2.1,0.5), limits=c(0,2.1))+
  xlab("Average dementia diagnosis delay (years)\n among those without cancer history")+ ylab("Observed IRR for dementia diagnosis \n(cancer history vs. no cancer history)")+
  geom_hline(yintercept=1, colour="black", lwd=1) +
  scale_color_manual(name="Interaction effect of cancer and dementia \non mortality rate (% of multiplicative)"
  ,values=brewer.pal(6,"Blues")[2:6], 
                     labels=c("40%", "50%", "60%", "70%", "80%"))+
  theme_bw()+
  theme(legend.position = c(0.2,0.2),
        axis.text.x = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_text(size=18, face="bold"), 
        axis.title.y = element_text(size=18, face="bold"), 
        legend.text = element_text(size=16)
  )
nodxbias

ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Figure5A.jpg", 
       device="jpg", plot=nodxbias, dpi="retina")



withdxbias<-ggplot(data=results[results$cancertype=="all" & results$Incr_contact==1 &results$Decr_demdx==0.8,])+
  geom_point(aes(x=unlist(Demdxwaittime), y=unlist(Demdx_IRR_PTcalc), 
                 group=factor(unlist(Mortmultfactor)), 
                 color=factor(unlist(Mortmultfactor))),
             size=2)+
  geom_line(aes(x=unlist(Demdxwaittime), y=unlist(Demdx_IRR_PTcalc), 
                group=factor(unlist(Mortmultfactor)), 
                color=factor(unlist(Mortmultfactor))), 
            size=2)+
  ylim(0.7,1.1)+
  scale_x_continuous(breaks=seq(0,2.1,0.5), limits=c(0,2.1))+
  xlab("Average dementia diagnosis delay (years)\n among those without cancer history")+ ylab("Observed IRR for dementia diagnosis \n(cancer history vs. no cancer history)")+
  geom_hline(yintercept=1, colour="black", lwd=1) +
  scale_color_manual(name="Interaction effect of cancer and dementia \non mortality rate (% of multiplicative)"
                     ,values=brewer.pal(6,"Blues")[2:6], 
                     labels=c("40%", "50%", "60%", "70%", "80%"))+
  theme_bw()+
  theme(legend.position = c(0.2,0.2),
        axis.text.x = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_text(size=18, face="bold"), 
        axis.title.y = element_text(size=18, face="bold"), 
        legend.text = element_text(size=16)
  )
withdxbias


ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Figure5B.jpg", 
       device="jpg", plot=withdxbias, dpi="retina")


#Sensitivity to dx bias parameters
results_dxbiasplot<-results[results$Mortmultfactor==0.5 & results$Demdxwaittime == 1.5,]

results_dxbiasplot$cancertype<-factor(results_dxbiasplot$cancertype, levels = c("breast", "prostate", "all", "lung"))

  dxbiassens<-ggplot(data=results_dxbiasplot)+
  geom_point(aes(x=unlist(dxbias), y=unlist(Demdx_IRR_PTcalc), 
                 group=factor(unlist(cancertype)), color=factor(unlist(cancertype))), size=2)+
  geom_line(aes(x=unlist(dxbias), y=unlist(Demdx_IRR_PTcalc), 
                group=factor(unlist(cancertype)), color=factor(unlist(cancertype))), size=2)+
  #ylim(0.7,1)+
  xlab("Relative rate of dementia diagnosis in those with vs. without cancer history")+ ylab("Observed IRR for dementia diagnosis \n(cancer history vs. no cancer history)")+
  geom_hline(yintercept=1, colour="black", lwd=1) +
    geom_vline(xintercept=1, colour="black", lwd=1) +
    theme_bw()+
  labs(color="Cancer type")+scale_color_discrete(name = "Cancer type", labels=c("Breast", "Prostate", "All", "Lung"))+
    scale_x_continuous(breaks=c(0.5, .7, .9, 1.0, 1.1, 1.3, 1.5))+
  theme(axis.text.x = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_text(size=18, face="bold"), 
        axis.title.y = element_text(size=18, face="bold"), 
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size=16),
        legend.title = element_text(size=20)
  )

  dxbiassens


ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Figure4.jpg", 
       device="jpg", plot=dxbiassens, dpi="retina")


#Framed as extra average delay


results_dxbiasplot$diff_delay<-((1/(results_dxbiasplot$dxbias*(1/as.numeric(results_dxbiasplot$Demdxwaittime))))-as.numeric(results_dxbiasplot$Demdxwaittime))*12
results_dxbiasplot$cancertype<-factor(results_dxbiasplot$cancertype, levels = c("breast", "prostate", "all", "lung"))


dxbiassens_v2<-ggplot(data=results_dxbiasplot)+
  geom_point(aes(x=unlist(diff_delay), y=unlist(Demdx_IRR_PTcalc), 
                 group=factor(unlist(cancertype)), color=factor(unlist(cancertype))), size=2)+
  geom_line(aes(x=unlist(diff_delay), y=unlist(Demdx_IRR_PTcalc), 
                group=factor(unlist(cancertype)), color=factor(unlist(cancertype))), size=2)+
  #ylim(0.7,1)+
  xlab("Difference in average dementia diagnosis delay between\n those with and without cancer (months)")+ ylab("Observed IRR for dementia diagnosis \n(cancer history vs. no cancer history)")+
  scale_x_continuous(breaks=seq(-5,15,5))+
  geom_hline(yintercept=1, colour="black", lwd=1) +
  geom_vline(xintercept=0, colour="black", lwd=1) +
  theme_bw()+
  labs(color="Cancer type")+scale_color_discrete(name = "Cancer type", labels=c("Breast", "Prostate", "All", "Lung"))+
  theme(axis.text.x = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_text(size=18, face="bold"), 
        axis.title.y = element_text(size=18, face="bold"), 
        legend.position = c(0.8, 0.7),
        legend.text = element_text(size=16),
        legend.title = element_text(size=20)
  )

dxbiassens_v2




ggsave("C:/Users/ehlarson/MHL Dropbox/Eleanor Hayes-Larson/UCLA/Eleanor_ERM/Cancer_AD sims/Github repo/Output/Dx bias figures/Figure4_v2.jpg", 
       device="jpg", plot=dxbiassens_v2, dpi="retina")
  


