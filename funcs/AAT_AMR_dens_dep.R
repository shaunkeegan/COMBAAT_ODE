## This is an ordinary differential equation model of African Animal 
## Trypanosomiasis (AAT) that incorporates the emergence, spread and loss of 
## antimicrobial resistance (AMR) between cattle, tsetse fly vectors and 
## wildlife. 

## VERSION 4 - May 2022
##
## This version includes compartments to account for the teneral phenomenon.


## Authors:   Shaun Keegan (shaun.keegan@glasgow.ac.uk)
##            Louise Matthews (louise.mattthews@glasgow.ac.uk)


## FORMAT: This file uses plain text descriptions of model parameters for user
##         accessibility. Mathematical model descriptions and corresponding 
##         parameter tables can be found at: 
##         http://github.com/shaunkeegan/AAT_AMR_main/model

## USAGE:  This file has been designed to be run and sourced from other files 
##         in the git repository, so that the model file is left untouched when 
##         exploring scenarios which are included at: 
##         http://github.com/shaunkeegan/AAT_AMR_main/scenarios

library(codetools)

AAT_AMR_dens_dep <- function(times, init, parms){
  
  # C - Cattle
  CS  <- init["CS"] # Susceptible
  CEs <- init["CEs"] # Exposed (drug sensitive strain)
  CEr <- init["CEr"] # Exposed (drug resistant strain)
  CIs <- init["CIs"] # Infected (drug sensitive strain)
  CIr <- init["CIr"] # Infected (drug resistant strain)
  CTs <- init["CTs"] # Treated (drug sensitive strain)
  CTr <- init["CTr"] # Treated (drug resistant strain)
  #CR  <- init["CR"] # Recovered
  
  # P - Prophylactically treated cattle
  PF <- init["PF"]  # Susceptible Fully protected
  PS  <- init["PS"]  # Susceptible
  PEs <- init["PEs"] # Exposed (drug sensitive strain)
  PEr <- init["PEr"] # Exposed (drug resistant strain)
  PIs <- init["PIs"] # Infected (drug sensitive strain)
  PIr <- init["PIr"] # Infected (drug resistant strain)
  PTs <- init["PTs"] # Treated (drug sensitive strain)
  PTr <- init["PTr"] # Treated (drug resistant strain)
  PPs <- init["PPs"] # Recovered
  PPr <- init["PPr"] # Recovered
  
  # W - Wildlife
  WS  <- init["WS"] # Susceptible
  WEs <- init["WEs"] # Exposed (drug sensitive strain)
  WEr <- init["WEr"] # Exposed (drug resistant strain)
  WIs <- init["WIs"] # Infected (drug sensitive strain)
  WIr <- init["WIr"] # Infected (drug resistant strain)
  #WR  <- init["WR"] # Recovered
  
  # V - Vectors
  VSt <- init["VSt"] # Susceptible teneral
  VSf <- init["VSf"] # Susceptible fed
  VEs <- init["VEs"] # Exposed (drug sensitive strain) 
  VEr <- init["VEr"] # Exposed (drug resistant strain)
  VIs <- init["VIs"] # Infected (drug sensitive strain)
  VIr <- init["VIr"] # Infected (drug resistant strain) 
  
  ## ----- Cattle
  birth.c          <- parms["birth.c"]
  biterate         <- parms["biterate"]
  prob.infection   <- parms["prob.infection"]
  gamma   <- parms["gamma.c"]
  death            <- parms["death.c"]
  sigma         <- parms["sigma.c"]
  treatment.q      <- parms["treatment.q"]
  treatment.p      <- parms["treatment.p"]
  sigma.st      <- parms["sigma.st"]
  emergence.p      <- parms["emergence.p"]  
  emergence.f      <- parms["emergence.f"]
  rec.adj          <- parms["rec.adj"]
  prop.prophylaxis <- parms["prop.prophylaxis"]
  fit.adj          <- parms["fit.adj"]
  waning           <- parms["waning"]
  new.prop         <- parms["new.prop"]    
  waning.f2s       <- parms["waning.f2s"]
  
  ## ----- Wildlife
  birth.w            <- parms["birth.w"]
  prob.infection.s.w <- parms["prob.infection.s.w"]
  prob.infection.r.w <- parms["prob.infection.r.w"]
  gamma.w   <- parms["gamma.w"]
  death.w            <- parms["death.w"]
  sigma.w         <- parms["sigma.w"]
  reversion          <- parms["reversion"]
  
  ## ----- Vectors
  K                <- parms["K"]
  feeding.rate     <-  parms["feeding.rate"]
  prob.infection.v <-  parms["prob.infection.v"]
  death.v <- parms["death.v"]
  birth.v <- parms["birth.v"]
  gamma.v <- parms["gamma.v"]
  ten2fed <- parms["ten2fed"]
  
  
  # Population total ----
  N <- CS + CEs + CEr + CIs + CIr + CTs + CTr +
    PF + PS + PEs + PEr + PIs + PIr + PTs + PTr +  PPs + PPr +
    WS + WEs + WEr + WIs + WIr 
  C <- CS + CEs + CEr + CIs + CIr + CTs + CTr 
  P <- PF + PS + PEs + PEr + PIs + PIr + PTs + PTr + PPs + PPr
  PC <- CS + CEs + CEr + CIs + CIr + CTs + CTr +
    PF + PS + PEs + PEr + PIs + PIr + PTs + PTr + PPs + PPr
  W <- WS + WEs + WEr + WIs + WIr 
  V <- VSt + VSf + VEs + VEr + VIs + VIr
  
  # Cattle ----
  # 
  # CS, CEs, CEr, CIs, CIr, CTs, CTr, CR
  
  dCS.dt <- birth.c * (1 - prop.prophylaxis) * PC +
    waning * PS -
    biterate * prob.infection * CS * VIs / N -  
    biterate * (prob.infection * fit.adj) * CS * VIr / N  + 
    sigma  * CIs + 
    sigma  * CIr + 
    sigma.st  * CTs +
#    sigma.st  * PPs +  #test addition
    (sigma * rec.adj)  * CTr - 
#    new.prop * CS - 
    death * CS 
  
  dCEs.dt <- biterate * prob.infection * CS * VIs / N - 
    gamma * CEs + 
    waning * PEs - 
    death * CEs 
  
  dCEr.dt <- biterate * (prob.infection * fit.adj) * CS * VIr / N - 
    gamma * CEr + 
    waning * PEr - 
    death * CEr
  
  dCIs.dt <- gamma * CEs - 
    treatment.q * CIs - 
    treatment.p * CIs - 
    sigma  * CIs + 
    waning * PIs + 
    waning * PPs - #LM moved from CTs equation
    death * CIs 
  
  dCIr.dt <- gamma * CEr - 
    treatment.q * CIr - 
    treatment.p * CIr - 
    sigma  * CIr + 
    waning * PIr +
    waning * PPr -  #LM moved from CTr equation 29/9/22
    death * CIr 
  
  dCTs.dt <- treatment.q * CIs - 
    sigma.st  * CTs - 
    emergence.f * CTs + 
    waning * PTs - 
#    waning * PPs - #LM moved up to CIs equation
    death * CTs
  
  dCTr.dt <- treatment.q * CIr - 
    (sigma * rec.adj) * CTr  + 
    emergence.f * CTs + 
    waning * PTr - 
#    waning *PPr - #LM moved up to CIr equation 29/9/22
    death * CTr
  
  
  
  # dCTs.dt <- 0
  # 
  # dCTr.dt <- 0
  # 
  # dCR.dt <- treatment.q * CIs + treatment.q * CIr + sigma  * CIs + sigma  * CIr + sigma  * CTs + 
  #   sigma  * CTr - resusceptible * CR - death * CR
  
  
  # Cattle with prophylaxis ----
  # 
  # PS, PEs, PEr, PIs, PIr, PTs, PTr, PR
  
  dPF.dt <- birth.c * (prop.prophylaxis) * PC - #new.prop * CS # Adding new prophylactically treated cattle
    biterate * (prob.infection * fit.adj * 1) * PF * VIr / N +   # Infection of resistant strain
    sigma.st  * PPs +                                     # sigma from treated (prophylactic) sensitive strain infection
    (sigma * rec.adj)  * PPr -                            # sigma from treated (prophylactic) resistant strain infection
    waning.f2s * PF -                                        # Waning prophylaxis from fully protected to partially protected
    death * PF                                               # Death of prophylactic susceptibles (fully protected)
  
  dPS.dt <-  waning.f2s * PF -                               # Waning of prophylactically treated cattle to semi protected
    biterate * prob.infection * PS * VIs / N -               # Infection of sensitive strain
    biterate * (prob.infection * fit.adj) * PS * VIr / N +   # Infection of resistant strain
    sigma  * PIs +                                        # sigma from sensitive strain infection
    sigma  * PIr +                                        # sigma from resistant strain infection
    sigma.st  * PTs +                                     # sigma from treated (fast acting) sensitive strain infection
    (sigma * rec.adj)  * PTr   -                          # sigma from treated (fast acting) resistant strain infection
    waning * PS -                                            # Waning of infection to non-prophylactic class
    death * PS                                               # Death of prophylactic susceptibles (partially protected)
  
  dPEs.dt <- biterate * prob.infection * PS * VIs / N -      # Infection of sensitive strain
    gamma * PEs -                                   # Movement from exposed to infectious
    emergence.p * PEs -
    waning *PEs -                                            # Waning of infection to non-prophylactic class
    death * PEs                                              # Death of prophylactic exposed (sensitive strain)
  
  dPEr.dt <- biterate * (prob.infection * fit.adj) * PS *VIr / N +    # Infection of resistant strain
    biterate * (prob.infection * fit.adj * 1) * PF * VIr / N -   # Infection of resistant strain
    gamma * PEr +                                   # Movement from exposed to infectious
    emergence.p * PEs -
    waning * PEr -                                           # Waning of infection to non-prophylactic class
    death * PEr                                              # Death of prophylactic exposed (resistant strain)
  
  dPIs.dt <- gamma * PEs -                                   # Movement from exposed to infectious
    treatment.q * PIs -                                      # Treatment with fast acting drug
    treatment.p * PIs -                                      # Treatment with prophylactic acting drug 
    sigma  * PIs -                                        # sigma from sensitive strain infection
    emergence.p * PIs -                                        # Emergence of AMR
    waning * PIs -                                           # Waning of infection to non-prophylactic class
    death * PIs                                              # Death of prophylactic infectious (sensitive strain)
  
  dPIr.dt <- gamma * PEr -                                   # Movement from exposed to infectious 
    treatment.q * PIr -                                      # Treatment with fast acting drug 
    treatment.p * PIr -                                      # Treatment with prophylactic acting drug 
    sigma  * PIr +                                        # sigma from resistant strain infection
    emergence.p * PIs -                                        # Emergence of AMR 
    waning * PIr -                                           # Waning of infection to non-prophylactic class
    death * PIr                                              # Death of prophylactic infectious (resistant strain)
  
  dPTs.dt <- treatment.q * PIs -                                      # Treatment with fast acting drug 
    sigma.st  * PTs -                                     # sigma from sensitive strain infection (fast acting treatment)
    emergence.p * PTs -    
    emergence.f * PTs -                                        # Emergence of AMR  
    waning * PTs -                                           # Waning of infection to non-prophylactic class 
    death * PTs                                              # Death of sensitive treated (fast acting)
  
  dPTr.dt <- treatment.q * PIr -                                      # Treatment with fast acting drug  
    (sigma * rec.adj)  * PTr +                            # Treatment with prophylactic acting drug  
    emergence.p * PTs +    
    emergence.f * PTs  -                                       # Emergence of AMR 
    waning * PTr -                                           # Waning of infection to non-prophylactic class  
    death * PTr                                              # Death of sensitive treated (prophylactic)  
  
  dPPs.dt <- treatment.p * PIs +                                      # Treatment with prophylactic acting drug  
    treatment.p * CIs -                                      # Treatment with prophylactic acting drug  
    emergence.p * PPs -
    sigma.st  * PPs -                                     # sigma from sensitive strain infection (prophylactic treatment) 
    waning * PPs -                                           # Waning of infection to non-prophylactic class 
    death* PPs                                               # Death of sensitive treated (prophylactic)  
  
  dPPr.dt <- treatment.p * PIr +                                      # Treatment with prophylactic acting drug   
    treatment.p * CIr +                                      # Treatment with prophylactic acting drug   
    emergence.p * PPs -
    (sigma * rec.adj)  * PPr -                            # sigma from resistant strain infection (prophylactic treatment)  
    waning * PPr -                                           # Waning of infection to non-prophylactic class 
    death* PPr                                               # Death of resistant treated (prophylactic) 
  
  
  # Wildlife ----
  # 
  # WS, WEs, WEr, WIs, WIr, WTs, WTr
  
  dWS.dt <- birth.w * W - biterate * prob.infection * WS * VIs / N - 
    biterate * (prob.infection * fit.adj) * WS * VIr / N  - 
    death.w * WS + sigma.w * WIs + sigma.w * WIr 
  
  dWEs.dt <- biterate * prob.infection * WS * VIs / N - gamma.w * 
    WEs - death.w * WEs 
  
  dWEr.dt <- biterate * (prob.infection * fit.adj) * WS * VIr / N - gamma.w * 
    WEr - death.w * WEr 
  
  dWIs.dt <- gamma.w * WEs - sigma.w * WIs - death.w * WIs + reversion * WIr
  
  dWIr.dt <- gamma.w * WEr - sigma.w * WIr - death.w * WIr - reversion * WIr
  
  # if (W < 1.0-8){
  #   dWS.dt <- 0
  #   dWEs.dt <- 0
  #   dWEr.dt <- 0
  #   dWIs.dt <- 0
  #   dWIr.dt <- 0
  # }
  # 
  
  # Tsetse ----
  # 
  # VS, VEs, VEr, VIs, VIr, 
  
  dVSt.dt <- birth.v * V  * (1 - V / K ) - 
    prob.infection.v * biterate * (CIs/N) * VSt -
    prob.infection.v * biterate * (CIr/N) * VSt -
    prob.infection.v * biterate * (CTs/N) * VSt -
    prob.infection.v * biterate * (CTr/N) * VSt -
    prob.infection.v * biterate * (PIs/N) * VSt -
    prob.infection.v * biterate * (PIr/N) * VSt -
    prob.infection.v * biterate * (PPs/N) * VSt -  #LM
    prob.infection.v * biterate * (PPr/N) * VSt -  #LM
    prob.infection.v * biterate * (PTs/N) * VSt -  #LM
    prob.infection.v * biterate * (PTr/N) * VSt -
    prob.infection.v * biterate * (WIs/N) * VSt -
    prob.infection.v * biterate * (WIr/N) * VSt -
#    prob.infection.v * biterate * (CIs + CIr + CTs +CTr + PIs + PIr + PPs + PPr + PTs + PTr + WIs + WIr) / N * VSt -
    ten2fed * VSt -
    death.v * VSt 
  
  dVSf.dt <- - 
    prob.infection.v * biterate * (CIs/N) * VSf - #LM missing minus at start of line
    prob.infection.v * biterate * (CIr/N) * VSf -
    prob.infection.v * biterate * (CTs/N) * VSf -
    prob.infection.v * biterate * (CTr/N) * VSf -
    prob.infection.v * biterate * (PIs/N) * VSf -
    prob.infection.v * biterate * (PIr/N) * VSf -
    prob.infection.v * biterate * (PPs/N) * VSf -  #LM addition
    prob.infection.v * biterate * (PPr/N) * VSf -  #LM addition
    prob.infection.v * biterate * (PTs/N) * VSf -  #LM addition
    prob.infection.v * biterate * (PTr/N) * VSf -
    prob.infection.v * biterate * (WIs/N) * VSf -
    prob.infection.v * biterate * (WIr/N) * VSf +
#  prob.infection.v * biterate * (CIs + CIr + CTs +CTr + PIs + PIr + PPs + PPr + PTs + PTr + WIs + WIr) / N * VSf +
    ten2fed * VSt -
    death.v * VSf
  
  dVEs.dt <-  + 
    prob.infection.v * biterate * (CIs/N) * VSt +
    prob.infection.v * biterate * (CTs/N) * VSt +
    prob.infection.v * biterate * (PIs/N) * VSt +
    prob.infection.v * biterate * (PPs/N) * VSt +  #LM addition
    prob.infection.v * biterate * (PTs/N) * VSt +
    prob.infection.v * biterate * (WIs/N) * VSt +
    prob.infection.v * biterate * (CIs/N) * VSf +
    prob.infection.v * biterate * (CTs/N) * VSf +
    prob.infection.v * biterate * (PIs/N) * VSf +
    prob.infection.v * biterate * (PPs/N) * VSf +  #LM addition
    prob.infection.v * biterate * (PTs/N) * VSf +
    prob.infection.v * biterate * (WIs/N) * VSf -
    #prob.infection.v * biterate * (CIs + 0*CIr + CTs + 0*CTr + PIs + 0*PIr + PPs + 0*PPr + PTs + 0*PTr + WIs + 0*WIr) / N * VSf +
    #prob.infection.v * biterate * (CIs + 0*CIr + CTs + 0*CTr + PIs + 0*PIr + PPs + 0*PPr + PTs + 0*PTr + WIs + 0*WIr) / N * VSt -
    gamma.v * VEs - death.v *VEs
  
  dVEr.dt <-  + 
    prob.infection.v * biterate * (CIr/N) * VSt +
    prob.infection.v * biterate * (CTr/N) * VSt +
    prob.infection.v * biterate * (PIr/N) * VSt +
    prob.infection.v * biterate * (PPr/N) * VSt +  #LM addition
    prob.infection.v * biterate * (PTr/N) * VSt +
    prob.infection.v * biterate * (WIr/N) * VSt +
    prob.infection.v * biterate * (CIr/N) * VSf +
    prob.infection.v * biterate * (CTr/N) * VSf +
    prob.infection.v * biterate * (PIr/N) * VSf +
    prob.infection.v * biterate * (PPr/N) * VSf +  #LM addition, corrected VSt to VSf 29/9/2022
    prob.infection.v * biterate * (PTr/N) * VSf +
    prob.infection.v * biterate * (WIr/N) * VSf -
#    prob.infection.v * biterate * (0*CIs + 1*CIr + 0*CTs + 1*CTr + 0*PIs + 1*PIr + 0*PPs + 1*PPr + 0*PTs + 1*PTr + 0*WIs + 1*WIr) / N * VSf +
#    prob.infection.v * biterate * (0*CIs + 1*CIr + 0*CTs + 1*CTr + 0*PIs + 1*PIr + 0*PPs + 1*PPr + 0*PTs + 1*PTr + 0*WIs + 1*WIr) / N * VSt -
    gamma.v * VEr - death.v *VEr
  
  dVIs.dt <- gamma.v * VEs - death.v * VIs
  
  dVIr.dt <- gamma.v * VEr - death.v * VIr
  
  # Model output ----
  #dX <- c(dCS.dt, dCEs.dt, dCEr.dt, dCIs.dt, dCIr.dt, dCTs.dt, dCTr.dt,  
  #        dPF.dt, dPS.dt, dPEs.dt, dPEr.dt, dPIs.dt, dPIr.dt, dPTs.dt, dPTr.dt, dPPs.dt, dPPr.dt,
  #        dWS.dt, dWEs.dt, dWEr.dt, dWIs.dt, dWIr.dt, 
  #        dVSt.dt, dVSf.dt, dVEs.dt, dVEr.dt, dVIs.dt, dVIr.dt)
  dX <- c(dCS.dt, dCEs.dt, 0.0, dCIs.dt, 0.0, dCTs.dt, 0.0,  
          dPF.dt, dPS.dt, dPEs.dt, 0.0, dPIs.dt, 0.0, dPTs.dt, 0.0, dPPs.dt, 0.0,
          dWS.dt, dWEs.dt, 0.0, dWIs.dt, 0.0, 
          dVSt.dt, dVSf.dt, dVEs.dt, 0.0, dVIs.dt, 0.0)
  list(dX)
  
  
}


findGlobals(fun = AAT_AMR_dens_dep, merge = FALSE)$variables
