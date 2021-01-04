#This script contains all of the functions for dementia incidence rates

# test<-rationalfit(AD_in_formod$Age[2:7], AD_in_formod$Rate_all[2:7], d1=4,d2=5)
# p1 <- test$p1; p2 <- test$p2
# ys <- polyval(p1,seq(0,40,1/12)) / polyval(p2,seq(0,40,1/12))
# plot(AD_in_formod$Age[2:7], AD_in_formod$Rate_all[2:7], type="l", col = "blue", xlim=c(0,40), ylim=c(0,0.3))
# points(seq(0,40,1/12), Re(ys), col="red")
# points(AD_in_formod$Age, AD_in_formod$Rate_all, col="blue", pch=15)


#Define reference dementia incidence rate as a function of time using ACT data
get.dem <- function(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                    calibration_data, Dem_inc){
  demratecorrection <- 
    if (S1.dem.rateratio == 1 | Sprev == 0) {as.matrix(rep(1, 6))
    } else{
      calibration_data[, paste0("demrate.adj.S_IRR", S1.dem.rateratio,
                                "_prev", Sprev, "_mortRR", S1.mort.rateratio)]
    }
  
  #I had to fill in NAs to put this in a dataframe so now get rid of them
  demratecorrection <- demratecorrection[!is.na(demratecorrection)]
  max_index = length(demratecorrection)
  
  if(time %% 5 == 0){
    index = time/5 + 1
  } else{
    index = ceiling(time/5)
  }
  
  if(index > max_index){index = max_index}
  
  AD_rate<- polyval(Dem_inc$p1,time) / polyval(Dem_inc$p2,time)
  
  return(Re(AD_rate)*demratecorrection[index])
} 

#------------------------------------------------------------------
# Defining dementia incidence rates based on time (i.e. age)
#------------------------------------------------------------------

#For those with survival characteristic, reference rate * protective IRR
get.dem.S1 <- function(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                       calibration_data, Dem_inc){
  return(get.dem(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                 calibration_data, Dem_inc)*S1.dem.rateratio)}

#For those without survival characteristic, dementia incidence rate is reference rate 
get.dem.S0 <- function(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                       calibration_data, Dem_inc){
  return(get.dem(time, S1.dem.rateratio, Sprev, S1.mort.rateratio, 
                 calibration_data, Dem_inc))}  


# #Unit test
# get.dem(40, 1,0,1,calibration_data,Dem_inc)
# get.dem.S1(6.5, 0.7,0,1,calibration_data,Dem_inc)
# get.dem.S0(6.5, 0.7,0,1,calibration_data,Dem_inc)
