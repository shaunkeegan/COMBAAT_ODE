##                                                                            ##
##                                                                            ##
##                    Explore Plot Functions // VERSION 2                     ##
##                                                                            ##
##                                                                            ##

## Required packages & files ----- 

library(ggplot2)
library(gghighlight)
library(cowplot)
library(dplyr)

source("funcs/plot_functions.R")

treatment_option <- "P"

if (treatment_option == "F"){
folder_name = "output/plots_fast_finer_insecticide/"
load("output/test_F.Rda")
}
if (treatment_option == "P"){
  folder_name = "output/plots_proph_finer_insecticide/"
  load("output/test_P.Rda")
}
if (treatment_option == "B"){
  folder_name = "output/plots_both_fit1/"
  load("output/test_B_fit1.Rda")
}

defaultW <- getOption("warn") 
options(warn = -1) 

## Load data and prep plot inputs -----

#load("output/test.Rda")
test1 <- test %>% signif(4)
p.ins <- unique(test1$prop.insecticide)
this_prop <- 0.05
ins2.dat <- test1 %>% filter(prop.insecticide == this_prop) 
this_prop <- 0.0
ins1.dat <- test1 %>% filter(prop.insecticide == this_prop) 
ins1.dat <- ins1.dat %>% filter(dose == 1)



plotlims.min <- c(0,0,0,0,0.6,0)
plotlims.max <- c(0.8,650,0.08,650,1, 65)
plotlims <- cbind(plotlims.min, plotlims.max)
names <- c("prevalence", "Incidence", "Risk", "No_trt_cat", "Prob_onward_tran", "R0_sen")
plotlims <- cbind(names, plotlims)
fitadj <- 0.95


## Run plots -----
pdf_width = 7
pdf_height = 6

## R0 sensitive (Y), number of wildlife (X), plot by treatment proportion
pdf(file= paste0(folder_name, "1_r0sen.pdf"), height = pdf_height, width = pdf_width)
plot_R0_Sen_dd(ins1.dat, c(0, 0.6,0.95), fitadj = fitadj)
dev.off()

pdf(file= paste0(folder_name, "2_resvsen.pdf"), height = pdf_height, width = pdf_width)
plot_res_v_sen_dd(ins1.dat, c(0,100,250), fitadj = fitadj)
dev.off()
## y options 
## R0_sen  R0_res R_eq_sen R_eq_res
## Prob_onward_tran Risk prevalence Incidence No_trt_cat
this_var <- "R0_sen"
y_label <- "R0 sensitive"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- as.numeric(axislims["plotlims.max"])
label <- "FALSE"
pdf(file= paste0(folder_name, "3_r0sen_tp.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, label = label, prev_threshold, ymin, ymax)
dev.off()

this_var <- "prevalence"
y_label <- "Prevalence"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- as.numeric(1)
label <- "FALSE"
pdf(file=paste0(folder_name, "4_prevsen_tp.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()

this_var <- "prevalence"
y_label <- "Prevalence"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- as.numeric(1)
label <- "FALSE"
pdf(file=paste0(folder_name, "4_prevsen_tp_v2.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()

this_var <- "Incidence"
y_label <- this_var
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- as.numeric(axislims["plotlims.max"])
label <- "FALSE"
pdf(file=paste0(folder_name, "5_incsen_tp.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()

this_var <- "No_trt_cat"
y_label <- "Number of Treated Cattle"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- as.numeric(axislims["plotlims.max"])
label <- "FALSE"
pdf(file= paste0(folder_name, "6_notrtsen_tp.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()


this_var <- "Prob_onward_tran"
y_label <- "Probability of Onward Transmission"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- as.numeric(axislims["plotlims.max"])
label <- "FALSE"
pdf(file=paste0(folder_name, "7_probontran_tp.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()


this_var <- "RiskA"
y_label <- "Relative risk of emergence"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9#1.1*max(ins1.dat$RiskA) #as.numeric(axislims["plotlims.max"])
label <- "FALSE"
#plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label, this_prop = this_prop)
pdf(file= paste0(folder_name, "8a_risk_Atp.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()

this_var <- "RiskE"
y_label <- "Relative risk of emergence and spread"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9#1.1*max(ins1.dat$RiskE) #as.numeric(axislims["plotlims.max"])
label <- "FALSE"
#plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label, this_prop = this_prop)
pdf(file= paste0(folder_name, "8_riskE_tp.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()




test1 <- test %>% signif(4)
p.ins <- unique(test1$prop.insecticide)
this_prop <- 0.0
ins1.dat <- test1 %>% filter(prop.insecticide == this_prop & treat_prop == 0.8) 


this_var <- "RiskE"
y_label <- "Relative risk of emergence"
prev_threshold <- 1
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9#1.1*max(ins1.dat$RiskA) #as.numeric(axislims["plotlims.max"])
label <- "FALSE"
#plot_PT.Y(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label, this_prop = this_prop)
pdf(file= paste0(folder_name, "dose_plot_tp0.8.pdf"), height = pdf_height, width = pdf_width)
plot_PT.Y.dose(ins1.dat, Wn = c(0,100,250), fitadj = fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label = label)
dev.off()
## Grid Plot -----

#Extra_plot
#source("Extra_plots.R")

##
## NEW MULTIPLOT //
## 
## 
##


K.10000 <- test1 %>% filter(K == 10000 & dose ==1, prop.insecticide < 0.5, prop.insecticide != 0.025) 
p.ins <- unique(K.10000$prop.insecticide)

this_var <- "RiskE"
y_label <- "Relative Risk of emergence"
threshold <- "prevalence"
prev_threshold <- 1
inc_threshold <- max(K.10000$Incidence)
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9*0.0001#1.1*max(ins1.dat$RiskA)*0.0001 #as.numeric(axislims["plotlims.max"])
col.var = "prop.insecticide"
shape.var = "W_st"
pdf(file=paste0(folder_name, "9a_10k_RiskE.pdf"))
plot_multiplot.insecticide(K.10000, Wn = c(0,100,250), p.ins,
                           col.var, shape.var,
                           fitadj = fitadj, my_theme, this_var, y_label, threshold, prev_threshold, inc_threshold, ymin, ymax , label = label)
dev.off()

#K.10000 <- test1 %>% filter(K == 10000 & dose ==1, prop.insecticide < 0.5, prop.insecticide != 0.025) 
this_var <- "RiskE"
y_label <- "Relative Risk of emergence"
prev_threshold <- 0.1
inc_threshold <- max(K.10000$Incidence)
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9*0.0001#1.1*max(ins1.dat$RiskA)*0.0001 #as.numeric(axislims["plotlims.max"])
col.var = "prop.insecticide"
shape.var = "W_st"
pdf(file=paste0(folder_name, "9b_10k_RiskE_pt0.1.pdf"))

plot_multiplot.insecticide(K.10000, Wn = c(0,100,250), p.ins,
                           col.var, shape.var,
                           fitadj = fitadj, my_theme, this_var, y_label, threshold, prev_threshold, inc_threshold, ymin, ymax , label = label)
dev.off()



#K.10000 <- test1 %>% filter(K == 10000 & dose == 1, prop.insecticide < 0.5, prop.insecticide != 0.025) 
this_var <- "RiskE"
y_label <- "Relative Risk of emergence"
prev_threshold <- 0.05
inc_threshold <- max(K.10000$Incidence)
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9*0.0001#1.1*max(ins1.dat$RiskA)*0.0001 #as.numeric(axislims["plotlims.max"])
col.var = "prop.insecticide"
shape.var = "W_st"
pdf(file=paste0(folder_name, "9c_10k_RiskE_pt0.05.pdf"))

plot_multiplot.insecticide(K.10000, Wn = c(0,100,250), p.ins,
                           col.var, shape.var,
                           fitadj = fitadj, my_theme, this_var, y_label, threshold, prev_threshold, inc_threshold, ymin, ymax , label = label)
dev.off()



#K.10000 <- test1 %>% filter(K == 10000& dose == 1, prop.insecticide < 0.5, prop.insecticide != 0.025) 
this_var <- "RiskE"
y_label <- "Relative Risk of emergence"
threshold <- "Incidence"
prev_threshold <- 1
inc_threshold <- max(K.10000$Incidence)
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9*0.0001#1.1*max(ins1.dat$RiskA)*0.0001 #as.numeric(axislims["plotlims.max"])
col.var = "prop.insecticide"
shape.var = "W_st"
pdf(file=paste0(folder_name, "9d_10k_RiskE_pt1.pdf"))

plot_multiplot.insecticide(K.10000, Wn = c(0,100,250), p.ins,
                           col.var, shape.var,
                           fitadj = fitadj, my_theme, this_var, y_label, threshold, prev_threshold, inc_threshold, ymin, ymax , label = label)
dev.off()

#K.10000 <- test1 %>% filter(K == 10000 & dose == 1, prop.insecticide < 0.5, prop.insecticide != 0.025) 
this_var <- "RiskE"
y_label <- "Relative Risk of emergence"
prev_threshold <- 0.1
inc_threshold <- 100#max(K.10000$Incidence)
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9*0.0001#1.1*max(ins1.dat$RiskA)*0.0001 #as.numeric(axislims["plotlims.max"])
col.var = "prop.insecticide"
shape.var = "W_st"
pdf(file=paste0(folder_name, "9e_10k_RiskE_pt0.1it100.pdf"))

plot_multiplot.insecticide(K.10000, Wn = c(0,100,250), p.ins,
                           col.var, shape.var,
                           fitadj = fitadj, my_theme, this_var, y_label, threshold, prev_threshold, inc_threshold, ymin, ymax , label = label)
dev.off()



#K.10000 <- test1 %>% filter(K == 10000 & dose == 1, prop.insecticide < 0.5, prop.insecticide != 0.025) 
this_var <- "RiskE"
y_label <- "Relative Risk of emergence"
prev_threshold <- 0.05
inc_threshold <- 10#max(K.10000$Incidence)
axislims <- plotlims[which (plotlims[,1]==this_var),]
ymin <- as.numeric(axislims["plotlims.min"])
ymax <- 9*0.0001#1.1*max(ins1.dat$RiskA)*0.0001 #as.numeric(axislims["plotlims.max"])
col.var = "prop.insecticide"
shape.var = "W_st"
pdf(file=paste0(folder_name, "9f_10k_RiskA_pt0.05it10.pdf"))

plot_multiplot.insecticide(K.10000, Wn = c(0,100,250), p.ins,
                           col.var, shape.var,
                           fitadj = fitadj, my_theme, this_var, y_label, threshold, prev_threshold, inc_threshold, ymin, ymax , label = label)
dev.off()






pdf(file = "My Plot.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 7) # The height of the plot in inches

# Step 2: Create the plot with R code
plot_R0_Sen_dd(ins1.dat, c(0,0.6,0.95), fitadj = fitadj)

# Step 3: Run dev.off() to create the file!
dev.off()
options(warn = defaultW)

source("Extra_plots_new.R")
