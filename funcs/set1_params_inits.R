## --------------------- Parms & Initis - set1
##
##
##
##





## ------------------------------------------------------ LOAD FUNCTIONS

library(codetools)


## Parameters & Initial Conditions ----

set1 <- function(output, birth.adj, fit.adj, K, prop_treat, prop.insecticide, NW, prop.prophylaxis, trt.type, dose.adj, emergence.adj){
  
  
  
  ## Cattle ----- 
  birth.c          <- 1 /(5*365)
  biterate         <- 0.8 / 4
  prob.infection   <- 0.46
  gamma.c   <- 1 / 15
  resusceptible    <- 10
  death.c            <- birth.c
  death.p          <- death.c
  sigma.c         <- 1 / 100 
  treatment       <- 1 * prop_treat * (sigma.c + death.c) / (1 - prop_treat)
  emergence       <- 0
  

  
  if(trt.type == "F"){
    treatment.q      <- treatment
    treatment.p      <- 0
    emergence.p      <- 0
    emergence.f      <- emergence * emergence.adj
  }else{
    if(trt.type == "P"){
    treatment.q      <- 0
    treatment.p      <- treatment
    emergence.p      <- emergence * emergence.adj
    emergence.f      <- 0
  }else{
    treatment.q      <- treatment
    treatment.p      <- treatment
    emergence.p      <- emergence * emergence.adj
    emergence.f      <- emergence * emergence.adj
  }}
  
  sigma.st      <- (1/3) * dose.adj + sigma.c * (1-dose.adj)    #* 250 #LM: adjusted so that R0 drops below 1 when 99% treated to reflect Hargrove
  rec.adj          <- 1
  waning           <- (1/30) / dose.adj
  waning.f2s       <- (1/60) / dose.adj
  new.prop         <- 0 
  
  
  NC <- 50    # Total cattle
  #figure out equilibrium in absence of infection
  #birth.c * prop.prophylaxis *NC - death.c * PF - waning.f2s*PF
  PF <- birth.c * prop.prophylaxis *NC / ( death.c + waning.f2s)
  #waning.f2s * PF - death.p * PS - waning * PS
  PS <- waning.f2s * PF / (death.p + waning)
  # birth.c * (1-prop.prophylaxis) * NC - death.c * CS + waning * PS
  CS <- (birth.c * (1-prop.prophylaxis) * NC + waning * PS)/death.c
  
  
  CIr <- 0    # Infected (drug resistant strain)
  CIs <- 1   # Infected (drug sensitive strain)
  #CS  <- NC * birth.c * (1 - prop.prophylaxis)/ death.c + waning  
  CS <- CS - CIs - CIr # Susceptible
  CEs <- 0    # Exposed (drug sensitive strain)
  CEr <- 0    # Exposed (drug resistant strain)
  
  CTs <- 0    # Treated (drug sensitive strain)
  CTr <- 0    # Treated (drug resistant strain)

  
  #PF  <- NC * birth.c * prop.prophylaxis/(death.c + waning.f2s)    # Susceptible   #LM: changed cattle to NC
  #PS  <- PF * waning.f2s/(death.p + waning)
  PEs <- 0    # Exposed (drug sensitive strain)
  PEr <- 0    # Exposed (drug resistant strain)
  PIs <- 0    # Infected (drug sensitive strain)
  PIr <- 0    # Infected (drug resistant strain)
  PTs <- 0    # Treated (drug sensitive strain)
  PTr <- 0    # Treated (drug resistant strain)
  PR  <- 0    # Recovered
  PPs <- 0
  PPr <- 0
  
  
  ## ----- Wildlife
  birth.w            <- 1 / 365
  prob.infection.s.w <- 0.46
  prob.infection.r.w <- 0.46
  gamma.w   <- 1 / 20
  resusceptible.w    <- 1 / 100
  death.w            <- birth.w
  sigma.w         <- sigma.c
  reversion          <- 0
  
  #NW <- 0
  WIs <- 0    # Infected (drug sensitive strain)
  WS  <- NW - WIs  # Susceptible
  WEs <- 0    # Exposed (drug sensitive strain)
  WEr <- 0    # Exposed (drug resistant strain)
  WIr <- 0    # Infected (drug resistant strain)

  
  
  ## -----  Vectors
  
  ten2fed <- 1/4
  qf <- 0.96 # Probability of surviving on a feeding day
  qn <- 0.98 # Probability of surviving on a non-feeding day
  feed.cyc <- 4 # Days between feeding
  feed.frequency     <-  0
  prob.infection.v <-  0.025
  incubation       <-  20
  prop.insecticide.actual <- prop.insecticide * NC/(NC + NW) # Proportion of insecticide adjusted for wildlife
  death.v <- -1 * log((1 - prop.insecticide.actual) * qf * qn ^ feed.cyc) / feed.cyc # Vector death rate
  birth.v = birth.adj *(-1) * log((1 - 0) * qf * qn ^ feed.cyc) / feed.cyc # Vector birth rate 
  equil_vector_pop <- max(0, K*(1 - death.v/birth.v)) # Vector equilibrium population
  gamma.v <- death.v * exp(-death.v * incubation) / (1 - exp(-death.v * incubation)) # Rate from E to I
  
  
  
  NV  <- equil_vector_pop  #equil_vector_pop
  VSt  <- NV  # Susceptible
  VSf <- 0
  VEs <- 0    # Exposed (drug sensitive strain)
  VEr <- 0    # Exposed (drug resistant strain)
  VIs <- 0    # Infected (drug sensitive strain)
  VIr <- 0    # Infected (drug resistant strain)
  
  
  
  
  
  ## ----- Parameters & initial conditions output
  
  params <- cbind(NC, NV, NW, birth.c, biterate, prob.infection, fit.adj, rec.adj, sigma.st, 
                  gamma.c, resusceptible, death.c, treatment.p, treatment.q, sigma.c, birth.v, 
                  death.v, feed.frequency, prob.infection.v, gamma.v, emergence.p, emergence.f,
                  reversion, K, birth.w, gamma.w, resusceptible.w, death.w, sigma.w, equil_vector_pop,
                  waning, waning.f2s, new.prop, ten2fed, prop.prophylaxis)
  names <- colnames(params)
  params <- as.vector(params)
  names(params) <- names
  
  inits <- cbind(CS, CEs, CEr, CIs, CIr, CTs, CTr, PF, PS, PEs, PEr, PIs,
                 PIr, PTs, PTr, PPs, PPr, WS, WEs, WEr, WIs, WIr, VSt, VSf, VEs, 
                 VEr, VIs, VIr)
  names <- colnames(inits)
  inits <- as.vector(inits)
  names(inits) <- names
  
  
  if (output == "P"){
    return(params)
  }else{
    return(inits)
  }
  
}


my_rootfun <- function(t, y, params) {
  return(c(y['CS'] - 20.0, y['CIs'] - 10.0))
}
my_rootfun2 <- function (t, y, params) {
  dstate <- unlist(AAT_AMR_dens_dep(t, y, params)) # rate of change vector
  condition1 <- (y['CIs'] - 1e-5)
  condition2 <- sum(abs(dstate)) - 1e-5
  return(c(condition1, condition2))
}



findGlobals(fun = set1, merge = FALSE)$variables

## Error Checks ----



#if(CIs + CIr > NC){cat(red("WARNING: Infected cattle (", CIs + CIr, ") is GREATER THAN total cattle(", NC, ")\n"))}
#if(WIs > NW){cat(red("WARNING: Infected wildlife (", WIs, ") is GREATER THAN total wildlife(", NW, ")\n"))}


