##                                                                            ##
##                                                                            ##
##                        Plot Functions // VERSION 2                         ##
##                                                                            ##
##                                                                            ##




## Plot theme -----

my_theme <- function () { 
  theme_grey(base_size=16) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(),
          axis.text.x = element_text(angle=45, hjust=1, size = 12),
          axis.text.y = element_text( size = 12),
          axis.title.x = element_text( size = 16),
          axis.title.y = element_text( size = 19),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

my_theme2 <- function () { 
  theme_grey(base_size=14) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(),
          axis.text.x = element_text(angle=45, hjust=1, size = 12),
          axis.title.x = element_text( size = 16),
          axis.title.y = element_text( size = 19),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

## Proportion Treated as X AXIS variable ----- 


plot_PT.Y <- function(df2, Wn, fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label, this_prop){
  
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  my_cols_vec <- ggplotColours(n=3)
  my_cols <- my_cols_vec[1:3]
  
  df2 <- df2 %>% filter(fit.adj == fitadj)
  names(df2)[names(df2) == this_var] <- "y"
  
  plot_this <- df2 %>% filter(W_st %in% Wn) %>% mutate(W_st = as.factor(W_st), K = as.factor(K))
  
  p <- ggplot(plot_this, aes(x = treat_prop, y = y, colour = W_st, group = interaction(W_st, K)) )  + my_theme() +
    geom_point(aes(shape = K), size = 3) + geom_line() +
    xlab("\n") + ylab("") +
    guides(col = guide_legend("Wildlife"), shape = guide_legend("Carrying \nCapacity (K)"))
  
  plots <- list()
  for (i in 1:3){
    plots[[i]] <- plot_this %>% filter(W_st == Wn[i]) %>% 
      ggplot(aes(x = treat_prop, y = y, group = interaction(W_st, K)) ) + my_theme() +
      geom_point(aes(shape = K), size = 3, colour = my_cols[i] ) + geom_line(colour = my_cols[i], size = 1.1  ) + ylim(ymin, ymax) +xlab("\n") + ylab("") +
      if (this_var != "prevalence"){gghighlight(prevalence < prev_threshold)}
  }
  
  legend <- get_legend(p + theme(legend.box.margin = margin(0, 0, 0, 12))) 
  
  if(label == "TRUE"){plot_grid(plots[[1]] + theme(legend.position="none") + ylab(y_label) +
                                  draw_figure_label(label = "Figure 1", fontface = "bold", position = "top.left"), 
                              plots[[2]] + theme(legend.position="none") + xlab("Proportion\nTreated"), 
                              plots[[3]] + theme(legend.position="none"), legend, ncol = 4)}else
                                {plot_grid(plots[[1]] + theme(legend.position="none") + ylab(y_label), 
                                 plots[[2]] + theme(legend.position="none") + xlab("Proportion\nTreated"), 
                                 plots[[3]] + theme(legend.position="none"), legend, ncol = 4)}
  
}



plot_multiplot.insecticide <- function(df2, Wn, p.ins, col.var, shape.var,
                                       fitadj, my_theme, this_var, y_label, threshold, prev_threshold, inc_threshold, ymin, ymax, label, this_prop){
  
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  my_cols_vec <- ggplotColours(n=5)
  my_cols <- my_cols_vec[1:5]
  
  df2 <- df2 %>% filter(fit.adj == fitadj)
  names(df2)[names(df2) == this_var] <- "y"
  
  plot_this <- df2 %>% filter(prop.insecticide %in% p.ins) %>% mutate(W_st = as.factor(W_st), prop.insecticide = as.factor(prop.insecticide))
  
  p <- ggplot(plot_this, aes(x = treat_prop, y = y*0.999, colour = prop.insecticide, group = interaction(W_st, prop.insecticide)) )  + my_theme() +
    geom_point(aes(shape = W_st), size = 2) + geom_line() +
    xlab("\n") + ylab("") +
    guides(col = guide_legend("Prop.insect"), shape = guide_legend("Wildlife"))
  
  plots <- list()
  for (i in 1:5){
    plots[[i]] <- plot_this %>% filter(prop.insecticide == p.ins[i]) %>% 
      ggplot(aes(x = treat_prop, y = y*0.999, group = interaction(W_st, prop.insecticide)) ) + my_theme() +
      geom_point(aes(shape = W_st), size = 2, colour = my_cols[i] ) + geom_line(colour = my_cols[i]  ) + ylim(ymin, ymax/0.0001) +xlab("\n") + ylab("")+ 
      if (threshold == "prevalence"){
        if (this_var != "prevalence"){gghighlight(prevalence < prev_threshold)}
      }else{gghighlight(Incidence < inc_threshold)}
      
  }
  
  legend <- get_legend(p + theme(legend.box.margin = margin(0, 0, 0, 12), legend.title=element_text(size=10)) )

  
  plot_grid(plots[[1]] + theme(legend.position="none") , 
            plots[[2]] + theme(legend.position="none"), 
            legend, 
            plots[[3]] + theme(legend.position="none")+ ylab(y_label), 
            plots[[4]] + theme(legend.position="none")+ xlab("Proportion\nTreated"),
            plots[[5]] + theme(legend.position="none"),ncol = 3)
  
}




## R0 sensitive plot -----

plot_R0_Sen_dd <- function(df2, trtprops, fitadj){
  
  my_par <- list(colour = "grey92")
  
  p <- ggplot(df2) + my_theme() +
    geom_point(aes(x = W_st, y = R0_sen, colour = as.factor(K)), size = 3) + 
    geom_line(aes(x = W_st, y = R0_sen, colour = as.factor(K))) +
    xlab("\n") +
    ylab("")+
    guides(col = guide_legend("Carrying \nCapacity (K)"))
  p
  p1 <- p + gghighlight(treat_prop == trtprops[1], fit.adj == fitadj,  unhighlighted_params = my_par) + ggtitle(paste0("TP = ",trtprops[1]))
  p2 <- p + gghighlight(treat_prop > trtprops[2]-0.0001, treat_prop < trtprops[2] + 0.0001, fit.adj == fitadj,   unhighlighted_params = my_par)+ ggtitle(paste0("TP = ",trtprops[2]))
  p3 <- p + gghighlight(treat_prop > trtprops[3]-0.0001, treat_prop < trtprops[3] +0.0001, fit.adj == fitadj,   unhighlighted_params = my_par) + ggtitle(paste0("TP = ",trtprops[3]))
  
  legend <- get_legend(p + theme(legend.box.margin = margin(0, 0, 0, 12), legend.title=element_text(size=15)) ) 
  plot_grid(p1 + theme(legend.position="none") + geom_line(aes(x = W_st, y = R0_sen, colour = as.factor(K),shape = as.factor(W_st))) + ylab("R0 sensitive"), 
            p2 + theme(legend.position="none") + geom_line(aes(x = W_st, y = R0_sen, colour = as.factor(K),shape = as.factor(W_st))) + xlab("Number of \nWildlife"), 
            p3 + theme(legend.position="none") + geom_line(aes(x = W_st, y = R0_sen, colour = as.factor(K),shape = as.factor(W_st))), legend, ncol = 4)
  
  
}


## R0 sensitive vs R0 resistant -----

plot_res_v_sen_dd <- function(df2, Wn, fitadj){
  
  #R0 plots: Resistant strain does better than sensitive apart from very low treatment rates
  my_par <- list(colour = "grey92")
  
  p1 <- df2 %>% filter(W_st %in% Wn) %>% filter(fit.adj %in% fitadj)  %>% ggplot() + my_theme() +
    geom_point(aes(x = treat_prop, y = R_eq_res/R_eq_sen, colour = as.factor(W_st), shape = as.factor(K)), size = 2) + 
    geom_line(aes(x = treat_prop, y = R_eq_res/R_eq_sen, colour = as.factor(W_st), shape = as.factor(K))) + 
    xlab("Treatment") + ylab("R resistant / R sensitive") + ylim(0,100)  + ggtitle(paste0("Fitness = ",fitadj)) +
    guides(col = guide_legend("Wildlife"), shape = guide_legend("Carrying \nCapacity (K)"))
  p1
  
  legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12))) 
  
  p2 <- df2 %>% filter(W_st %in% Wn) %>% filter(fit.adj %in% fitadj)  %>%ggplot() + my_theme() +
    geom_point(aes(x = treat_prop, y = R_eq_res/R_eq_sen, colour = as.factor(W_st), shape = as.factor(K)), size = 2) + 
    geom_line(aes(x = treat_prop, y = R_eq_res/R_eq_sen, colour = as.factor(W_st), shape = as.factor(K))) + 
    geom_abline(intercept = 1, slope = 0, linetype = 2) + 
    xlab("Treatment") + ylab("R resistant / R sensitive") + ylim(0,2)  + ggtitle(paste0("Fitness = ",fitadj)) +
    guides(col = guide_legend("Wildlife"), shape = guide_legend("Vectors"))
  p2 #+ gghighlight(R_eq_res/R_eq_sen < 1)
  plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), legend, ncol = 3)
}



plot_PT.Y.dose <- function(df2, Wn, fitadj, my_theme, my_par, this_var, y_label, prev_threshold, ymin, ymax, label, this_prop){
  
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  my_cols_vec <- ggplotColours(n=3)
  my_cols <- my_cols_vec[1:3]
  
  df2 <- df2 %>% filter(fit.adj == fitadj)
  names(df2)[names(df2) == this_var] <- "y"
  
  plot_this <- df2 %>% filter(W_st %in% Wn) %>% mutate(W_st = as.factor(W_st), K = as.factor(K))
  
  p <- ggplot(plot_this, aes(x = dose, y = y, colour = W_st, group = interaction(W_st, K)) )  + my_theme() +
    geom_point(aes(shape = K), size = 3) + geom_line() +
    xlab("\n") + ylab("") +
    guides(col = guide_legend("Wildlife"), shape = guide_legend("Carrying \nCapacity (K)"))
  
  plots <- list()
  for (i in 1:3){
    plots[[i]] <- plot_this %>% filter(W_st == Wn[i]) %>% 
      ggplot(aes(x = dose, y = y, group = interaction(W_st, K)) ) + my_theme() +
      geom_point(aes(shape = K), size = 3, colour = my_cols[i] ) + geom_line(colour = my_cols[i]  ) + ylim(ymin, ymax) +xlab("\n") + ylab("") +
      if (this_var != "prevalence"){gghighlight(prevalence < prev_threshold)}
  }
  
  legend <- get_legend(p + theme(legend.box.margin = margin(0, 0, 0, 12))) 
  
  if(label == "TRUE"){plot_grid(plots[[1]] + theme(legend.position="none") + ylab(y_label) +
                                  draw_figure_label(label = "Figure 1", fontface = "bold", position = "top.left"), 
                                plots[[2]] + theme(legend.position="none") + xlab("Dose\n"), 
                                plots[[3]] + theme(legend.position="none"), legend, ncol = 4)}else
                                {plot_grid(plots[[1]] + theme(legend.position="none") + ylab(y_label), 
                                           plots[[2]] + theme(legend.position="none") + xlab("Dose\n"), 
                                           plots[[3]] + theme(legend.position="none"), legend, ncol = 4)}
  
}



