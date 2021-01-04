#This script contains all of the functions for dementia incidence rates

#Define reference dementia incidence rate as a function of time using ACT data
get.demdx <- function(time, Demdx_inc){
  # demdxratecorrection <- 
  #   if (S1.dem.rateratio == 1 | Sprev == 0) {as.matrix(rep(1, 6))
  #   } else{
  #     calibration_data[, paste0("demrate.adj.S_IRR", S1.dem.rateratio,
  #                               "_prev", Sprev, "_mortRR", S1.mort.rateratio)]
  #   }
  # 
  # #I had to fill in NAs to put this in a dataframe so now get rid of them
  # demdxratecorrection <- demratecorrection[!is.na(demratecorrection)]
  max_index = 8
  # 
  if(time %% 5 == 0){
    index = time/5 + 1
  } else{
    index = ceiling(time/5)
  }
  
  if(index > max_index){index = max_index}
  
  return(Demdx_inc[index])#*demdxratecorrection[index])
} 



#------------------------------------------------------------------
# Defining dementia incidence rates based on time (i.e. age)
#------------------------------------------------------------------

#For those with cancer, reference rate * increase in contact with medical 
#          system * reduction in dx IRR
get.demdx.C1 <- function(time, Incr_contact, Decr_demdx, 
                         Demdx_inc){
  return(get.demdx(time, Demdx_inc)*Incr_contact*Decr_demdx)}

#For those without cancer, dementia incidence rate is reference rate 
get.demdx.C0 <- function(time, Demdx_inc){
  return(get.demdx(time, Demdx_inc))}  