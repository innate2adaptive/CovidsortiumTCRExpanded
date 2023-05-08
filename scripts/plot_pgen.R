# HCW alpha pgen

library(ggplot2)
library(plyr)
library(data.table)

dfa<-read.csv("data/output_data/exp_A_wide5_pgen.csv")
dfa0<-dfa[(dfa$pgen>0) & (dfa$control == "False"),]
# remove und because we are not interested
dfa1<-dfa0[dfa0$max_timepoint_class!="und",]

dfa_nonexp<-fread("data/output_data/HCW_A_controls_emerson_pGen_calculated.csv.gz")
dfa_nonexp$max_timepoint_class<-"non_expanded"
dfa_nonexp1<-dfa_nonexp[(dfa_nonexp$control==FALSE) & (dfa_nonexp$pgen>0),]

DFA<-rbind(dfa1[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")], 
      dfa_nonexp1[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")])
DFA<-DFA[!duplicated(DFA[,c("v_call", "junction_aa", "max_timepoint_class")]),] # to make it comparable to emerson

DFA$pgen<-as.numeric(as.character(DFA$pgen))
DFA$pgen_log10<-round(log10(DFA$pgen))
DFA$pgen_log10<-as.numeric(as.character(DFA$pgen_log10))
DFA$max_timepoint_class<-factor(DFA$max_timepoint_class, 
                                levels=c("non_expanded", "early", "late"))


dfa.new<-ddply(DFA,.(max_timepoint_class),summarise,
               prop=prop.table(table(pgen_log10)),
               pgen_log10=names(table(pgen_log10)))
dfa.new$pgen_log10<-as.numeric(as.character(dfa.new$pgen_log10))
dfa.new$prop<-as.numeric(as.character(dfa.new$prop))
dfa.new$max_timepoint_class<-factor(dfa.new$max_timepoint_class, 
                                levels=c("non_expanded", "early", "late"))
sum(dfa.new$prop)

p0<-ggplot(dfa.new,aes(x = pgen_log10, y=prop, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_col(position=position_dodge(preserve="single", width=0.8), width=0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gainsboro')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gray')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_HCWa.svg")
print(p0)
dev.off()

p01<-ggplot(dfa.new,aes(x = pgen_log10, y=prop, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_line(size=0.8) + geom_point(size=3) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gainsboro')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gray')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_HCWa_lines.svg")
print(p01)
dev.off()

p02<-ggplot(dfa.new[dfa.new$max_timepoint_class!="non_expanded",],aes(x = pgen_log10, y=prop, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_col(position=position_dodge(preserve="single", width=0.8), width=0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3")) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3")) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_HCWa_no-controls.svg")
print(p02)
dev.off()


dfb<-read.csv("data/output_data/exp_B_wide5_pgen.csv")
dfb<-dfb[(dfb$pgen>0) & (dfb$control == "False"),]
dfb1<-dfb[dfb$max_timepoint_class!="und",]

dfb_nonexp<-fread("data/output_data/HCW_B_controls_emerson_pGen_calculated.csv.gz")
dfb_nonexp$max_timepoint_class<-"non_expanded"
dfb_nonexp<-dfb_nonexp[(dfb_nonexp$control==FALSE) & (dfb_nonexp$pgen>0),]

DFB<-rbind(dfb1[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")], 
           dfb_nonexp[,c("v_call", "junction_aa", "ID", "pgen", "max_timepoint_class")])
DFB<-DFB[!duplicated(DFB[,c("v_call", "junction_aa", "max_timepoint_class")]),] # to make it comparable to emerson

DFB$pgen<-as.numeric(as.character(DFB$pgen))
DFB$pgen_log10<-round(log10(DFB$pgen))
DFB$pgen_log10<-as.numeric(as.character(DFB$pgen_log10))
DFB$max_timepoint_class<-factor(DFB$max_timepoint_class, 
                                levels=c("non_expanded", "early", "late"))

dfb.new<-ddply(DFB,.(max_timepoint_class),summarise,
               prop=prop.table(table(pgen_log10)),
               pgen_log10=names(table(pgen_log10)))
dfb.new$pgen_log10<-as.numeric(as.character(dfb.new$pgen_log10))
dfb.new$prop<-as.numeric(as.character(dfb.new$prop))
dfb.new$max_timepoint_class<-factor(dfb.new$max_timepoint_class, 
                                    levels=c("non_expanded", "early", "late"))

p1<-ggplot(dfb.new,aes(x = pgen_log10, y=prop, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_col(position=position_dodge(preserve="single", width=0.8), width=0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gainsboro')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gray')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_HCWb.svg")
print(p1)
dev.off()

p11<-ggplot(dfb.new,aes(x = pgen_log10, y=prop, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_line(size=0.8) + geom_point(size=3) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gainsboro')) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3", non_expanded='gray')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_HCWb_lines.svg")
print(p11)
dev.off()

p12<-ggplot(dfb.new[dfb.new$max_timepoint_class!="non_expanded",],aes(x = pgen_log10, y=prop, fill=max_timepoint_class, color=max_timepoint_class))+
  geom_col(position=position_dodge(preserve="single", width=0.8), width=0.5) +
  scale_fill_manual(values = c(early = "darkorange2", late = "cyan3")) +
  scale_color_manual(values = c(early = "darkorange2", late = "cyan3")) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_HCWb_no-controls.svg")
print(p12)
dev.off()

lcmv<-read.csv("data/output_data/LCMV_pgen.csv")
lcmv<-lcmv[lcmv$pgen>0,]
lcmv<-lcmv[lcmv$condition != "PBS",] # we are not interested in these
lcmv<-lcmv[!duplicated(lcmv[,c("junction_aa","day")]), ]# to make it comparable to precursor analysis
lcmv_c<-read.csv("data/output_data/LCMV_controls_pgen.csv")
lcmv_c$condition<-"controls"
lcmv_c<-lcmv_c[lcmv_c$pgen>0,]
lcmv_c<-lcmv_c[!duplicated(lcmv_c[,c("junction_aa","day")]), ]# to make it comparable to precursor analysis
lcmv1<-rbind(lcmv[,c("junction_aa", "day", "pgen")], 
             lcmv_c[,c("junction_aa", "day", "pgen")])
lcmv1$pgen<-as.numeric(as.character(lcmv1$pgen))
lcmv1$pgen_log10<-round(log10(lcmv1$pgen))
lcmv1$pgen_log10<-as.numeric(as.character(lcmv1$pgen_log10))

lcmv1$day<-factor(lcmv1$day, 
                       levels=c("controls","day8", "day40"))

lcmv.new<-ddply(lcmv1,.(day),summarise,
             prop=prop.table(table(pgen_log10)),
             pgen_log10=names(table(pgen_log10)))
lcmv.new$pgen_log10<-as.numeric(as.character(lcmv.new$pgen_log10))
lcmv.new$prop<-as.numeric(as.character(lcmv.new$prop))
lcmv.new$day<-factor(lcmv.new$day, 
                                     levels=c("controls", "day8", "day40"))

p2<-ggplot(lcmv.new,aes(x = pgen_log10, y=prop, fill=day, color=day))+
  geom_col(position=position_dodge(preserve="single", width=0.8), width=0.5) +
  scale_fill_manual(values = c(day8 = "darkorange2", day40 = "cyan3", controls='gainsboro')) +
  scale_color_manual(values = c(day8 = "darkorange2", day40 = "cyan3", controls='gray')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_LCMV.svg")
print(p2)
dev.off()

p21<-ggplot(lcmv.new,aes(x = pgen_log10, y=prop, fill=day, color=day))+
  geom_line(size=0.8) + geom_point(size=3) +
  scale_fill_manual(values = c(day8 = "darkorange2", day40 = "cyan3", controls='gainsboro')) +
  scale_color_manual(values = c(day8 = "darkorange2", day40 = "cyan3", controls='gray')) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_LCMV_lines.svg")
print(p21)
dev.off()

p22<-ggplot(lcmv.new[lcmv.new$day!="controls",],aes(x = pgen_log10, y=prop, fill=day, color=day))+
  geom_col(position=position_dodge(preserve="single", width=0.8), width=0.5) +
  scale_fill_manual(values = c(day8 = "darkorange2", day40 = "cyan3")) +
  scale_color_manual(values = c(day8 = "darkorange2", day40 = "cyan3")) +
  xlim(-23,-4) +
  labs(x = "log10(pgen)", y = "proportion of V/CDR3s", fill = "", color = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/pgen_LCMV_no-controls.svg")
print(p22)
dev.off()

print(p0)
print(p1)
print(p2)

print(p01)
print(p11)
print(p21)

