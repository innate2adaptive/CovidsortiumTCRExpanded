## this script plots boxplots from bulk analysis
## legacy bulk_analysis_plot_expanded.R from Martina's git

library(ggplot2)
library(reshape)

alpha_ss<-readRDS("data/output_data/exp_A_bulk_analysis_subsampled.RData")
beta_ss<-readRDS("data/output_data/exp_B_bulk_analysis_subsampled.RData")

alpha_ss$week_group<-unlist(lapply(as.numeric(as.character(alpha_ss$week)), function(x){
  ifelse(x < 0, "before", ifelse(x < 14, "during", "after"))
}))

alpha_ss$week<-factor(alpha_ss$week,levels=c(-3, -2, -1, 0, 1, 2, 3, 4, 14))
alpha_ss$week_group<-factor(alpha_ss$week_group,levels=c("before", "during", "after"))

mycols<-c("forestgreen", "darkseagreen3", "darkseagreen2", "red", "blue4", "cornflowerblue", "cyan3", "cyan", "darkgoldenrod")
names(mycols)<-c(-3, -2, -1, 0, 1, 2, 3, 4, 14)

alpha_ss$mygroup<-paste(alpha_ss$id, alpha_ss$week, sep = "-")

library(lemon)

control_names<-c("FALSE" = "PCR+", "TRUE" = "PCR-")

for (metric in c("richness_dcb", "shannon_dcb")){
  
  nc <- alpha_ss[alpha_ss$control == FALSE,]
  
  # print("non controls")
  print(metric)
  print("week -1, 0")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 0,][metric])))
  print(p$p.value)
  print("week -1, 1")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 1,][metric])))
  print(p$p.value)
  print("week -1, 2")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 2,][metric])))
  print(p$p.value)
  print("week -1, 3")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 3,][metric])))
  print(p$p.value)
  print("week -1, 4")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 4,][metric])))
  print(p$p.value)
  print("week -1, 14")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 14,][metric])))
  print(p$p.value)
  
  
  p1<- ggplot(alpha_ss, aes(x = week, y = as.numeric(.data[[metric]]))) +
    geom_line(aes(group = id), col = "grey", alpha = 0.5) +
    geom_point() + 
    scale_colour_manual(values = mycols) + 
    geom_boxplot(fill = NA)+
    scale_x_discrete(labels = c("-1", "0", "1", "2", "3", "4", "12-14")) +
    facet_rep_grid(rows = vars(control), labeller = as_labeller(control_names)) +
    theme_classic() + labs(y = strsplit(metric, "_")[[1]][1]) + theme(aspect.ratio = 0.8)
  svg(file=paste('output_figures/Fig1b_expanded_alpha_subsample_', metric, ".svg", sep = ""), width=5,height=5)
  print(p1)
  dev.off()
}

beta_ss$week_group<-unlist(lapply(as.numeric(as.character(beta_ss$week)), function(x){
  ifelse(x < 0, "before", ifelse(x < 14, "during", "after"))
}))
beta_ss$week<-factor(beta_ss$week,levels=c(-3, -2, -1, 0, 1, 2, 3, 4, 14))
beta_ss$week_group<-factor(beta_ss$week_group,levels=c("before", "during", "after"))


beta_ss$mygroup<-paste(beta_ss$id, beta_ss$week, sep = "-")

for (metric in c("richness_dcb", "shannon_dcb")){
  
  c <- beta_ss[beta_ss$control == TRUE,]
  nc <- beta_ss[beta_ss$control == FALSE,]
  
  # print("non controls")
  print(metric)
  print("week -1, 0")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 0,][metric])))
  print(p$p.value)
  print("week -1, 1")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 1,][metric])))
  print(p$p.value)
  print("week -1, 2")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 2,][metric])))
  print(p$p.value)
  print("week -1, 3")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 3,][metric])))
  print(p$p.value)
  print("week -1, 4")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 4,][metric])))
  print(p$p.value)
  print("week -1, 14")
  p<-wilcox.test(as.numeric(unlist(nc[nc$week == -1,][metric])), as.numeric(unlist(nc[nc$week == 14,][metric])))
  print(p$p.value)
  p1<- ggplot(beta_ss, aes(x = week, y = as.numeric(.data[[metric]]))) +
    geom_line(aes(group = id), col = "grey", alpha = 0.5) +
    geom_point() + 
    scale_colour_manual(values = mycols) + 
    geom_boxplot(fill = NA)+
    scale_x_discrete(labels = c("-1", "0", "1", "2", "3", "4", "12-14")) +
    facet_rep_grid(rows = vars(control), labeller = as_labeller(control_names)) +
    theme_classic() + labs(y = strsplit(metric, "_")[[1]][1]) + theme(aspect.ratio = 0.8)
  svg(file=paste('output_figures/Fig1b_expanded_beta_subsample_', metric, ".svg", sep = ""), width=5,height=5)
  print(p1)
  dev.off()
}

