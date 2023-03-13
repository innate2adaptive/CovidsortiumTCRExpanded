# On samples that have all timepoints, annotate early vs late

library(reshape)
library(ggplot2)

load("data/output_data/exp_AB_wide3.RData")

max_col <- unlist(max.col(exp_AB_wide3[,c('proportion_0', 'proportion_1', 'proportion_2', 'proportion_3', 'proportion_4')]))
exp_AB_wide3['max_timepoint']<-c('proportion_0', 'proportion_1', 'proportion_2', 'proportion_3', 'proportion_4')[max_col]
exp_AB_wide3['max_timepoint_class']<-unlist(lapply(exp_AB_wide3$max_timepoint, function(x){
                                            ifelse(x %in% c('proportion_0', 'proportion_1'), "early",
                                            ifelse(x %in% c('proportion_4', 'proportion_4'), "late",
                                            "und"))}))
hist(exp_AB_wide3[exp_AB_wide3$max_timepoint_class == 'late',]$proportion_0, 
     breaks = seq(from=0, to=410, by=10))
hist(exp_AB_wide3[exp_AB_wide3$max_timepoint_class == 'late',]$proportion_1, 
     breaks = seq(from=0, to=700, by=10))

exp_AB_wide3[exp_AB_wide3$max_timepoint_class == 'late' & (exp_AB_wide3$proportion_0>50 | exp_AB_wide3$proportion_1>50),]$max_timepoint_class <- 'und'
exp_AB_wide3[exp_AB_wide3$max_timepoint_class == 'early' & exp_AB_wide3$proportion_0<50 & exp_AB_wide3$proportion_1<50,]$max_timepoint_class <- 'und'

exp_AB_wide3[exp_AB_wide3$control == T, ]$max_timepoint_class <- 'und'

exp_melt<-melt(data.frame(exp_AB_wide3[!is.na(exp_AB_wide3$max_timepoint),c('proportion_0', 'proportion_1', 'proportion_2', 'proportion_3', 'proportion_4', 'proportion_14', 'max_timepoint_class', 'junction', 'ID')]),
              id = c('junction', 'ID', 'max_timepoint_class'))

exp_melt1<-exp_melt[exp_melt$max_timepoint_class != 'und',]

p<-ggplot(exp_melt1) +
  geom_jitter(aes(x=variable, y = value, col = max_timepoint_class), alpha = 0.2, position = position_dodge(width = 0.1), size = 3) +
  geom_line(aes(x=variable, y = value, col = max_timepoint_class, group = junction), alpha = 0.1) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", und='grey')) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(labels = c(0,1,2,3,4, 14)) +
  labs(x = "week since first PCR+", y = "frequency in sample", color = "max expansion") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        title=element_text(size=12),
        axis.text.x = element_text(vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)+
  facet_grid(cols=vars(max_timepoint_class))

svg("output_figures/HCW_early_late_definition.svg")
print(p)
dev.off()

save(exp_AB_wide3,file=  "data/output_data/exp_AB_wide4.RData")
