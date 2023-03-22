# HCW alpha pgen

library(ggplot2)
library(dplyr)

dfa<-read.csv("data/output_data/exp_A_wide5_pgen.csv")
dfa<-dfa[dfa$pgen>0,]
dfa<-dfa[!duplicated(dfa[,c("v_call", "junction_aa")]),]
dfa$pgen<-as.numeric(as.character(dfa$pgen))
dfa$pgen_log10<-round(log10(dfa$pgen))
dfa$pgen_log10<-as.numeric(as.character(dfa$pgen_log10))
dfa$max_timepoint_class<-factor(dfa$max_timepoint_class, 
                                levels=c("und", "early", "late"))

dfa.new<-ddply(dfa,.(max_timepoint_class),summarise,
               prop=prop.table(table(pgen_log10)),
               pgen_log10=names(table(pgen_log10)))
dfa.new$pgen_log10<-as.numeric(as.character(dfa.new$pgen_log10))
dfa.new$max_timepoint_class<-factor(dfa.new$max_timepoint_class, 
                                levels=c("und", "early", "late"))

ggplot(dfa.new,aes(x = pgen_log10, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_bar(aes(y=prop), stat="identity", position="identity", alpha = 0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", und='grey')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", und='grey')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)


dfb<-read.csv("data/output_data/exp_B_wide5_pgen.csv")
dfb<-dfb[dfb$pgen>0,]
dfb<-dfb[!duplicated(dfb[,c("v_call", "junction_aa")]),]
dfb$pgen<-as.numeric(as.character(dfb$pgen))
dfb$pgen_log10<-round(log10(dfb$pgen))
dfb$pgen_log10<-as.numeric(as.character(dfb$pgen_log10))
dfb$max_timepoint_class<-factor(dfb$max_timepoint_class, 
                                levels=c("und", "early", "late"))

dfb.new<-ddply(dfb,.(max_timepoint_class),summarise,
               prop=prop.table(table(pgen_log10)),
               pgen_log10=names(table(pgen_log10)))
dfb.new$pgen_log10<-as.numeric(as.character(dfb.new$pgen_log10))
dfb.new$max_timepoint_class<-factor(dfb.new$max_timepoint_class, 
                                    levels=c("und", "early", "late"))

ggplot(dfb.new,aes(x = pgen_log10, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_bar(aes(y=prop), stat="identity", position="identity", alpha = 0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", und='grey')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", und='grey')) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  xlim(-23,-4) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

lcmv<-read.csv("data/output_data/LCMV_pgen.csv")
lcmv<-lcmv[lcmv$pgen>0,]
lcmv$pgen<-as.numeric(as.character(lcmv$pgen))
lcmv$pgen_log10<-round(log10(lcmv$pgen))
lcmv$pgen_log10<-as.numeric(as.character(lcmv$pgen_log10))
lcmv$condition<-factor(lcmv$condition, 
                                levels=c("PBS", "LCMV8", "LCMV40"))

lcmv.new<-ddply(lcmv,.(condition),summarise,
               prop=prop.table(table(pgen_log10)),
               pgen_log10=names(table(pgen_log10)))
lcmv.new$pgen_log10<-as.numeric(as.character(lcmv.new$pgen_log10))
lcmv.new$condition<-factor(lcmv.new$condition, 
                                     levels=c("PBS", "LCMV8", "LCMV40"))

ggplot(lcmv.new,aes(x = pgen_log10, fill=condition, color=condition))+
  geom_bar(aes(y=prop), stat="identity", position="identity", alpha = 0.5) +
  scale_fill_manual(values = c(LCMV8 = "darkorange2", LCMV40 = "cyan3", PBS='grey'))+
  scale_color_manual(values = c(LCMV8 = "darkorange2", LCMV40 = "cyan3", PBS='grey')) +
  xlim(-23,-3) +
  labs(x = "log10(pgen)", y = "proportion of CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)
