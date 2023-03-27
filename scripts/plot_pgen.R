# HCW alpha pgen

library(ggplot2)
library(plyr)
library(data.table)

dfa<-read.csv("data/output_data/exp_A_wide5_pgen.csv")
dfa<-dfa[(dfa$pgen>0) & (dfa$control == "False"),]
dfa<-dfa[!duplicated(dfa[,c("v_call", "junction_aa", "max_timepoint_class")]),] # to make it comparable to emerson
# remove und because we are not interested
dfa1<-dfa[dfa$max_timepoint_class!="und",]

dfa_nonexp<-fread("data/output_data/HCW_A_controls_emerson_pGen_calculated.csv.gz")
dfa_nonexp$max_timepoint_class<-"non.expanded"

DFA<-rbind(dfa1[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")], 
      dfa_nonexp[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")])

DFA$pgen<-as.numeric(as.character(DFA$pgen))
DFA$pgen_log10<-round(log10(DFA$pgen))
DFA$pgen_log10<-as.numeric(as.character(DFA$pgen_log10))
DFA$max_timepoint_class<-factor(DFA$max_timepoint_class, 
                                levels=c("non.expanded", "early", "late"))


dfa.new<-ddply(DFA,.(max_timepoint_class),summarise,
               prop=prop.table(table(pgen_log10)),
               pgen_log10=names(table(pgen_log10)))
dfa.new$pgen_log10<-as.numeric(as.character(dfa.new$pgen_log10))
dfa.new$max_timepoint_class<-factor(dfa.new$max_timepoint_class, 
                                levels=c("non.expanded", "early", "late"))

p<-ggplot(dfa.new,aes(x = pgen_log10, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_bar(aes(y=prop), stat="identity", position="identity", alpha = 0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", non.expanded='lightgrey')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", non.expanded='lightgrey')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_HCWa.svg")
print(p)
dev.off()

dfb<-read.csv("data/output_data/exp_B_wide5_pgen.csv")
dfb<-dfb[(dfb$pgen>0) & (dfb$control == "False"),]
dfb<-dfb[!duplicated(dfb[,c("v_call", "junction_aa", "max_timepoint_class")]),] # to make it comparable to emerson
dfb1<-dfb[dfb$max_timepoint_class!="und",]

dfb_nonexp<-fread("data/output_data/HCW_B_controls_emerson_pGen_calculated.csv.gz")
dfb_nonexp$max_timepoint_class<-"non.expanded"

DFB<-rbind(dfb1[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")], 
           dfb_nonexp[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")])


DFB$pgen<-as.numeric(as.character(DFB$pgen))
DFB$pgen_log10<-round(log10(DFB$pgen))
DFB$pgen_log10<-as.numeric(as.character(DFB$pgen_log10))
DFB$max_timepoint_class<-factor(DFB$max_timepoint_class, 
                                levels=c("und", "early", "late"))

dfb.new<-ddply(DFB,.(max_timepoint_class),summarise,
               prop=prop.table(table(pgen_log10)),
               pgen_log10=names(table(pgen_log10)))
dfb.new$pgen_log10<-as.numeric(as.character(dfb.new$pgen_log10))
dfb.new$max_timepoint_class<-factor(dfb.new$max_timepoint_class, 
                                    levels=c("non.expanded", "early", "late"))

p1<-ggplot(dfb.new,aes(x = pgen_log10, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_bar(aes(y=prop), stat="identity", position="identity", alpha = 0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", non.expanded='lightgrey')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", non.expanded='lightgrey')) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  xlim(-23,-4) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)
svg("output_figures/pgen_HCWb.svg")
print(p1)
dev.off()

lcmv<-read.csv("data/output_data/LCMV_pgen.csv")
lcmv<-lcmv[lcmv$pgen>0,]
lcmv<-lcmv[!duplicated(lcmv[,c("junction_aa","condition")]), ]# to make it comparable to precursor analysis
lcmv<-lcmv[lcmv$condition != "PBS",] # we are not interested in these
lcmv_c<-read.csv("data/output_data/LCMV_controls_pgen.csv")
lcmv_c$condition<-"controls"
lcmv1<-rbind(lcmv[,c("junction_aa", "condition", "pgen")], 
             lcmv_c[,c("junction_aa", "condition", "pgen")])
lcmv1$pgen<-as.numeric(as.character(lcmv1$pgen))
lcmv1$pgen_log10<-round(log10(lcmv1$pgen))
lcmv1$pgen_log10<-as.numeric(as.character(lcmv1$pgen_log10))

lcmv1$condition<-factor(lcmv1$condition, 
                       levels=c("controls","LCMV8", "LCMV40"))

lcmv.new<-ddply(lcmv1,.(condition),summarise,
             prop=prop.table(table(pgen_log10)),
             pgen_log10=names(table(pgen_log10)))
lcmv.new$pgen_log10<-as.numeric(as.character(lcmv.new$pgen_log10))
lcmv.new$condition<-factor(lcmv.new$condition, 
                                     levels=c("controls", "LCMV8", "LCMV40"))

p2<-ggplot(lcmv.new,aes(x = pgen_log10, fill=condition, color=condition))+
  geom_bar(aes(y=prop), stat="identity", position="identity", alpha = 0.5) +
  scale_fill_manual(values = c(LCMV8 = "darkorange2", LCMV40 = "cyan3", controls='lightgrey'))+
  scale_color_manual(values = c(LCMV8 = "darkorange2", LCMV40 = "cyan3",controls='lightgrey')) +
  xlim(-23,-3) +
  labs(x = "log10(pgen)", y = "proportion of CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)
svg("output_figures/pgen_LCMV.svg")
print(p2)
dev.off()

print(p)
print(p1)
print(p2)

