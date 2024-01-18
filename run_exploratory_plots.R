## --------------------- EXECUTE: exploratory_plots
##
##
##
##





## ------------------------------------------------------ LOAD FUNCTIONS

library(deSolve)
library(tictoc)
library(progress)
library(ggplot2)
library(dplyr)
library(gghighlight)
library(cowplot)
library(crayon)

source("funcs/AAT_AMR_dens_dep.R")
source("funcs/r0.R")
source("funcs/set1_params_inits.R")
source("funcs/qual_check.R")
source("funcs/quick_plot.R")
source("funcs/loop_model_run.R")

loops <- TRUE


if (loops == TRUE) {
  treatment_type = "P"  #F this means quick treatment
  N_wl <- c(0, 100, 250)
  treat.prop.vecA <- seq(0, 0.9, by = 0.2)      #full from 0-1
  treat.prop.vecB <- seq(0.91, 0.99, by = 0.02)
  treat.prop.vec <- c(treat.prop.vecA, treat.prop.vecB)
  K.vec <- c(10000, 6000, 2000, 1000, 500)
  fit.adj.vec <- c(0.95)
  birth.adj.vec <- c(2)
  prop.insecticide.vec <- c(0.0,0.025, 0.05, 0.10,0.15, 0.2, 0.3, 0.5, 0.8,0.99)
  prop.prophylaxis <- 0.0
  dose.adj.vec <- 1#c(1, 0.9, 0.7, 0.4, 0.2)
  emergence.adj.vec <- 1
} else {
  treatment_type = "F" 
  N_wl <- 250
  treat.prop.vec <- 0.99
  #equil.vector_pop <- c(3000)
  fit.adj.vec <- c(0.6)
  K.vec <- 2000
  prop.insecticide.vec <- 0.0
  birth.adj.vec <- 2
  prop.prophylaxis <- 0.5
  dose.adj.vec <- 0.1
  emergence.adj.vec <- 1
}


df2 <- model_run(treatment_type, N_wl, treat.prop.vec, K.vec, fit.adj.vec, birth.adj.vec,
         prop.insecticide.vec, prop.prophylaxis, dose.adj.vec, emergence.adj.vec)





test <- df2
time <- format(Sys.time(), "%a %b %d %X %Y")
save(test, file = paste0("output_spk/test_", treatment_type, ".Rda"))
#save(test,file ="output/test.Rda")

quick_plot(out)


