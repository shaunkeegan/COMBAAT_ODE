## --------------------- R0
##
##
##
##




## ------------------------------------------------------ Sensitive Strain

r0_sen_calc <- function(params, inits){
  
  Nc <- inits["CS"]  
  Np <- inits["PS"]
  Nw <- inits["WS"]
  Nv <- inits["VSt"]
  Nh <- as.numeric(Nc + Np + Nw)
  
  biterate <- params["biterate"]
  
  treatment.p <- params["treatment.p"]
  treatment.q <- params["treatment.q"]
  waning <- params["waning"]
  
  gamma.c <- params["gamma.c"]
  death.c <- params["death.c"]
  sigma.c <- params["sigma.c"]
  
  gamma.p <- params["gamma.c"]
  death.p <- params["death.c"]
  sigma.p <- params["sigma.c"]
  
  gamma.w <- params["gamma.w"]
  death.w <- params["death.w"]
  sigma.w <- params["sigma.w"]
  
  gamma.v <- params["gamma.v"]
  death.v <- params["death.v"]
  
  RCV <- biterate * Nc / Nh * gamma.c / (gamma.c + death.c)
  RCV <- as.numeric(RCV)
  
  RPV <- biterate * Np / Nh * gamma.p / (gamma.p + death.p)
  RPV <- as.numeric(RPV)
  
  RWV <- biterate * Nw / Nh * gamma.w / (gamma.w + death.w)
  RWV <- as.numeric(RWV)
  
  # Probability of I -> Tp
  p1c <- treatment.p/ (treatment.p + treatment.q + sigma.c + death.c) 
  
  # Probability of Tp -> I
  p2c <- waning/ (treatment.p + treatment.q + sigma.c + death.c) 
  
  
  RVC <- biterate * Nv / Nh * gamma.v / (gamma.v + death.v) *
    ( 1/(treatment.p + treatment.q + sigma.c + death.c) * 1/(1-p1c *p2c) +
        treatment.q/(treatment.p + treatment.q + sigma.c + death.c) * 1/(death.c + sigma.c) +
        treatment.p/(treatment.p + treatment.q + sigma.c + death.c) * 1/ (death.c + sigma.c + waning) *
        1/(1-p1c *p2c))
  RVC <- as.numeric(RVC)
  
  # Probability of I -> Tp
  p1p <- treatment.p/ (treatment.p + treatment.q + sigma.p + death.p) 
  
  # Probability of Tp -> I
  p2p <- waning/ (treatment.p + treatment.q + sigma.p + death.p) 
  
  
  RVP <- biterate * Nv / Nh * gamma.v / (gamma.v + death.v) * 
    ( 1/(treatment.p + treatment.q + sigma.p + death.p) * 1/(1-p1c *p2c) +
        treatment.q/(treatment.p + treatment.q + sigma.p + death.p) * 1/(death.p + sigma.p) +
        treatment.p/(treatment.p + treatment.q + sigma.p + death.p) * 1/ (death.p + sigma.p + waning) *
        1/(1-p1c *p2c))
  RVP <- as.numeric(RVP)

  
  RVW <- biterate * Nv / Nh * 1 / (sigma.w + death.w) * gamma.v / (gamma.v + death.v)
  RVW <- as.numeric(RVW)
  

  R0 <- RCV * RVC + RPV * RVP + RWV * RVW
  
  names <- c("R0", "RCV", "RPV", "RWV", "RVC", "RVP", "RVW")
  output <- c(R0, RCV, RPV, RWV, RVC, RVP, RVW)
  names(output) <- names
  
  return(output)
  
}


r0_sen_calc_LM2 <- function(params, Nc, Np, Nw, Nv){
  
  Nh <- params["NC"] + params["NW"]  
  
  biterate <- params["biterate"]
  prob.infection <- params["prob.infection"]
  prob.infection.v <- params["prob.infection.v"]
  
  
  treatment.p <- params["treatment.p"]
  treatment.q <- params["treatment.q"]
  waning <- params["waning"]
  
  gamma.c <- params["gamma.c"]
  death.c <- params["death.c"]
  sigma.c <- params["sigma.c"]
  sigma.st <- params["sigma.st"]
  
  gamma.p <- params["gamma.c"]
  death.p <- params["death.c"]
  sigma.p <- params["sigma.c"]
  
  gamma.w <- params["gamma.w"]
  death.w <- params["death.w"]
  sigma.w <- params["sigma.w"]
  
  gamma.v <- params["gamma.v"]
  death.v <- params["death.v"]
  
  RCV <- biterate * prob.infection * Nc / Nh * gamma.c / (gamma.c + death.c) * 1 / (death.v)
  RCV <- as.numeric(RCV)
  
  RPV <- biterate * prob.infection * Np / Nh * gamma.p / (gamma.p + death.p) * 1 / (death.v)
  RPV <- as.numeric(RPV)
  
  RWV <- biterate * prob.infection * Nw / Nh * gamma.w / (gamma.w + death.w) * 1 / (death.v)
  RWV <- as.numeric(RWV)
  
  # Probability of I -> Tp
  p1c <- treatment.p/ (treatment.p + treatment.q + sigma.c + death.c)  
  
  # Probability of Tp -> I
  #p2c <- waning/ (treatment.p + treatment.q + sigma.st + death.c) 
  p2c <- waning/ (waning + sigma.st + death.c) #LM corrected
  
  
  RVC <- biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * 
    ( 1/(treatment.p + treatment.q + sigma.c + death.c) * 1/(1-p1c *p2c)  +
        treatment.q/(treatment.p + treatment.q + sigma.c + death.c) * 1/(death.c + sigma.st) * 1/(1-p1c *p2c) +
        treatment.p/(treatment.p + treatment.q + sigma.c + death.c) * 1/ (death.c + sigma.st + waning) *
        1/(1-p1c *p2c))
  RVC <- as.numeric(RVC)
  
  
  RVP <- biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * 
    (1/(treatment.p + treatment.q + sigma.p + death.p + waning)) +              #contribution from PIs
    
    (waning /(treatment.p + treatment.q + sigma.p + death.p + waning)) * RVC  + #contribution from waning back to CIS
    
    biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * (
      treatment.q/(treatment.p + treatment.q + sigma.p + death.p + waning) * ( 1/(death.p + sigma.st + waning) + #contribution from PTs
                                                                                 waning/(death.p + sigma.st + waning) * 1/(sigma.st + death.c)))  + # waning from PTs back to CTs 
    
    biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * (
      treatment.p /(treatment.p + treatment.q + sigma.p + death.p + waning)) * 1/(sigma.st + death.p + waning) + #contrib from PPs
    
    treatment.p /(treatment.p + treatment.q + sigma.p + death.p + waning) *
    waning /(sigma.st + death.p + waning) * RVC #contribution from PPs waning back to CIs
  
  
  RVP <- as.numeric(RVP)
  
  
  RVW <- biterate * prob.infection.v * Nv / Nh * 1 / (sigma.w + death.w) * gamma.v / (gamma.v + death.v)
  RVW <- as.numeric(RVW)
  
  
  R0 <- RCV * RVC + RPV * RVP + RWV * RVW
  
  names <- c("R0", "RCV", "RPV", "RWV", "RVC", "RVP", "RVW")
  output <- c(R0, RCV, RPV, RWV, RVC, RVP, RVW)
  names(output) <- names
  output
  
  return(output)
  
}


r0_calc_sen_or_res <- function(params, Nc, Np, Nw, Nv, sen){
  
  Nh <- params["NC"] + params["NW"]  
  
  biterate <- params["biterate"]
  prob.infection <- params["prob.infection"]
  prob.infection.v <- params["prob.infection.v"]
  fit.adj <- params["fit.adj"]
  
  treatment.p <- params["treatment.p"]
  treatment.q <- params["treatment.q"]
  waning <- params["waning"]
  
  gamma.c <- params["gamma.c"]
  death.c <- params["death.c"]
  sigma.c <- params["sigma.c"]
  sigma.st <- params["sigma.st"]
  
  gamma.p <- params["gamma.c"]
  death.p <- params["death.c"]
  sigma.p <- params["sigma.c"]
  
  gamma.w <- params["gamma.w"]
  death.w <- params["death.w"]
  sigma.w <- params["sigma.w"]
  
  gamma.v <- params["gamma.v"]
  death.v <- params["death.v"]
  
  
  if (sen == "yes"){sigma.treated <- sigma.st}
  if (sen == "no"){sigma.treated <- sigma.c}
  if (sen == "no"){prob.infection <- prob.infection* fit.adj}
  
  RCV <- biterate * prob.infection * Nc / Nh * gamma.c / (gamma.c + death.c) * 1 / (death.v)
  RCV <- as.numeric(RCV)
  
  RPV <- biterate * prob.infection * Np / Nh * gamma.p / (gamma.p + death.p) * 1 / (death.v)
  RPV <- as.numeric(RPV)
  
  RWV <- biterate * prob.infection * Nw / Nh * gamma.w / (gamma.w + death.w) * 1 / (death.v)
  RWV <- as.numeric(RWV)

  
  # Probability of I -> Tp
  p1c <- treatment.p/ (treatment.p + treatment.q + sigma.c + death.c)  
  
  # Probability of Tp -> I
  #p2c <- waning/ (treatment.p + treatment.q + sigma.st + death.c) 
  p2c <- waning/ (waning + sigma.treated + death.c) #LM corrected
  
  
  RVC <- biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * 
    ( 1/(treatment.p + treatment.q + sigma.c + death.c) * 1/(1-p1c *p2c)  +
        treatment.q/(treatment.p + treatment.q + sigma.c + death.c) * 1/(death.c + sigma.treated) * 1/(1-p1c *p2c) +
        treatment.p/(treatment.p + treatment.q + sigma.c + death.c) * 1/ (death.c + sigma.treated + waning) *
        1/(1-p1c *p2c))
  RVC <- as.numeric(RVC)
  
  
  RVP <- biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * 
    (1/(treatment.p + treatment.q + sigma.p + death.p + waning)) +              #contribution from PIs
    
    (waning /(treatment.p + treatment.q + sigma.p + death.p + waning)) * RVC  + #contribution from waning back to CIS
    
    biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * (
      treatment.q/(treatment.p + treatment.q + sigma.p + death.p + waning) * ( 1/(death.p + sigma.treated + waning) + #contribution from PTs
                                                                                 waning/(death.p + sigma.treated + waning) * 1/(sigma.treated + death.c)))  + # waning from PTs back to CTs 
    
    biterate * prob.infection.v * Nv / Nh * gamma.v / (gamma.v + death.v) * (
      treatment.p /(treatment.p + treatment.q + sigma.p + death.p + waning)) * 1/(sigma.treated + death.p + waning) + #contrib from PPs
    
    treatment.p /(treatment.p + treatment.q + sigma.p + death.p + waning) *
    waning /(sigma.treated + death.p + waning) * RVC #contribution from PPs waning back to CIs
  
  
  RVP <- as.numeric(RVP)
  
  
  RVW <- biterate * prob.infection.v * Nv / Nh * 1 / (sigma.w + death.w) * gamma.v / (gamma.v + death.v)
  RVW <- as.numeric(RVW)
  
  
  R0 <- RCV * RVC + RPV * RVP + RWV * RVW
  
  names <- c("R0", "RCV", "RPV", "RWV", "RVC", "RVP", "RVW")
  output <- c(R0, RCV, RPV, RWV, RVC, RVP, RVW)
  names(output) <- names
  output
  
  return(output)
  
}

