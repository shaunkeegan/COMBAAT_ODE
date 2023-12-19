## --------------------- Quality Check - qual_check
##
##
##
##





## ------------------------------------------------------ Check Inputs

library(crayon)
library(dplyr) #added by LM


qual_check_no0 <- function(input){
  
  test <- c(1:length(input))
  test[] <- NA
  
  for (i in 1: length(input)){
    if(input[i] < 0){cat(red("WARNING: ", names(input[i]), "  LESS THAN 0\n"))}
  }
  
  for(j in 1:length(input)){
    if (input[j] >= 0){test[j] <- 1}
  }

    if(sum(test[!is.na(test)]) == length(input)){cat(green("CHECK OK\n"))}
  
}



## ------------------------------------------------------ Check R = 1


equilibrium_R <- function(R0input, ODEinput){
  
  RCV <- R0input["RCV"] 
  RPV <- R0input["RPV"] 
  RWV <- R0input["RWV"] 
  RVC <- R0input["RVC"] 
  RVP <- R0input["RVP"] 
  RVW <- R0input["RVW"] 

  
  

  ODEinput2 <- tail(ODEinput,1)
  
  fc <- ODEinput2$CS / (ODEinput2$CS + ODEinput2$CEs + ODEinput2$CEr + 
                          ODEinput2$CIs + ODEinput2$CIr + ODEinput2$CTs +
                          ODEinput2$CTr)
  
  
  fp <- (ODEinput2$PF +  ODEinput2$PS) / (ODEinput2$PF + ODEinput2$PS + 
                                            ODEinput2$PEs + ODEinput2$PEr + 
                                            ODEinput2$PIs + ODEinput2$PIr + 
                                            ODEinput2$PTs + ODEinput2$PTr +
                                            ODEinput2$PPr + ODEinput2$PPs)
  fw <- ODEinput2$WS / (ODEinput2$WS + ODEinput2$WEs + ODEinput2$WEr + 
                          ODEinput2$WIs + ODEinput2$WIr)
  fv <- (ODEinput2$VSt + ODEinput2$VSf) / (ODEinput2$VSt + ODEinput2$VSf + 
                                             ODEinput2$VEs + ODEinput2$VEr + 
                                             ODEinput2$VIs + ODEinput2$VIr)
  
  
  
  
  return(as.numeric(RCV * fc * RVC * fv + RPV * fp * RVP * fv + RWV * fw * RVW * fv))
  
}


## ------------------------------------------------------ Add total cols to out


totals <- function(ODEinput){
  
  
  cattle.total <- ODEinput$CS + ODEinput$CEs + ODEinput$CEr + ODEinput$CIs +
                    ODEinput$CIr + ODEinput$CTs + ODEinput$CTr
  prophylactic.total <-ODEinput$PF + ODEinput$PS + ODEinput$PEs + ODEinput$PEr +
                        ODEinput$PIs + ODEinput$PIr + ODEinput$PTs +
                         ODEinput$PTr + ODEinput$PPr + ODEinput$PPs
  all.cows <- cattle.total+prophylactic.total
  wildlife.total <- ODEinput$WS + ODEinput$WEs + ODEinput$WEr +
                     ODEinput$WIs + ODEinput$WIr
  vector.total <- ODEinput$VSt + ODEinput$VSf + ODEinput$VEs + ODEinput$VEr +
                    ODEinput$VIs + ODEinput$VIr


  temp1 <- cbind(cattle.total, prophylactic.total, all.cows, wildlife.total, vector.total)
  ODEoutput <- cbind(ODEinput,temp1)
  
  return(ODEoutput)
  
}


totals_LM <- function(ODEinput){
  
  ODEoutput <- ODEinput %>% mutate(cattle.total = rowSums(select(., starts_with("C"))),
                                   prophylactic.total = rowSums(select(., starts_with("P"))),
                                   vector.total = rowSums(select(., starts_with("V"))),
                                   wildlife.total = rowSums(select(., starts_with("W"))),
                                   all.cows = cattle.total + prophylactic.total)
  
  return(ODEoutput)
  
}
