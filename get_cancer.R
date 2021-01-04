#------------------------------------------------------------------
# Defining cancer incidence rates based on time (i.e. age)
#------------------------------------------------------------------
#Cancer incidence rates with artificial final rates
# SEER_in_formod<-SEER_inc_obs[,c("time","Rate_all")]
# SEER_in_formod$time<-(SEER_in_formod$time-2.5)
# SEER_in_formod
# add_rows<-data.frame(time=c(27.5, 32.5, 37.5), Rate_all=SEER_in_formod$Rate_all[SEER_in_formod$time==22.5])
# SEER_in_formod<-rbind(SEER_in_formod, add_rows)
# SEER_in_formod
# test2<-rationalfit(SEER_in_formod$time, SEER_in_formod$Rate_all, d1=2,d2=3)
# p1 <- test2$p1; p2 <- test2$p2
# ys <- polyval(p1,seq(0,40,1/12)) / polyval(p2,seq(0,40,1/12))

#PLot to check if shape looks good. 
# plot(SEER_in_formod$time, SEER_in_formod$Rate_all, type="l", col = "blue", xlim=c(0,40), ylim=c(0,0.05))
# points(seq(0,40,1/12), Re(ys), col="red")
# points(SEER_in_formod$time, SEER_in_formod$Rate_all, col="blue", pch=15)


#Define cancer incidence rate as a function of time using SEER data
get.cancer <- function(time, Cancer_inc){
  
  
  SEER_rate<- polyval(Cancer_inc$p1,time) / polyval(Cancer_inc$p2,time)

  return(Re(SEER_rate))
}

#unit test
#get.cancer(16,Cancer_inc)
