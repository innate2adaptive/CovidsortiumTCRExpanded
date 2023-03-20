# load and compute emerson sharing

library(ggplot2)
library(reshape)
library(dplyr)

sharing_ems<-readRDS("data/downloaded_data/Emerson_sharing_aa.rds")
load("data/output_data/exp_AB_wide3.RData")

exp_b<-exp_AB_wide3[(exp_AB_wide3$chain == "beta"),] # emerson set only has beta
# exp_b<-exp_b[!duplicated(exp_b[,c("junction_aa")]),]

rm(exp_AB_wide3)

# merge with sharing - adds a column which says how many emerson individuals the cdr3 is found in
sharing_exp<-merge(exp_b, sharing_ems, by.x = "junction_aa", by.y = "aminoAcid", all.x = TRUE) # add sharing info
sharing_exp[is.na(sharing_exp$sharing_level),]$sharing_level<-0
dim(sharing_exp[sharing_exp$sharing_level >= 2,]) # 2903 as in paper

rm(exp_b)

load("data/control_beta_1.RData")
ctrl<-data.frame(unlist(control_b))
colnames(ctrl)<-c("junction_aa")
# ctrl<-data.frame(ctrl[!duplicated(ctrl[,c("junction_aa")]),])
# colnames(ctrl)<-c("junction_aa")

# merge with sharing - adds a column which says how many emerson individuals the cdr3 is found in
sharing_ctrl_B<-merge(ctrl, sharing_ems, by.x = "junction_aa", by.y = "aminoAcid", all.x = TRUE) # add sharing info
sharing_ctrl_B[is.na(sharing_ctrl_B$sharing_level),]$sharing_level<-0
dim(sharing_ctrl_B[sharing_ctrl_B$sharing_level >= 2,]) # 

sharing_ctrl_B$set<-"ctrl"
sharing_exp$set<-"exp"

sharing<-rbind(sharing_ctrl_B[c("junction_aa", "sharing_level", "set")], sharing_exp[c("junction_aa", "sharing_level", "set")])

counts<-data.frame(table(sharing$sharing_level, sharing$set))
counts$class <- unlist(lapply(as.numeric(as.character(counts$Var1)), function(x) {ifelse(x <= 3, "0-3",
                                                                                         ifelse((x>3) & (x <= 10), "4-10", 
                                                                                                ifelse((x>10) & (x <= 30), "11-30",
                                                                                                       ifelse((x>30) & (x <= 100), "31-100", 
                                                                                                              ifelse((x>100) & (x <= 300), "101-300", "301+")))))}))

sum_ctrl <- sum(as.numeric(as.character(counts[counts$Var2 == "ctrl",]$Freq)))
sum_exp <- sum(as.numeric(as.character(counts[counts$Var2 == "exp",]$Freq)))

counts$prop<-0

counts[counts$Var2 == "ctrl",]$prop <- counts[counts$Var2 == "ctrl",]$Freq/sum_ctrl
counts[counts$Var2 == "exp",]$prop <- counts[counts$Var2 == "exp",]$Freq/sum_exp

counts$class <- factor(counts$class, levels = c("0-3", "4-10", "11-30", "31-100", "101-300", ">301"))

counts_agg<-aggregate(counts$prop, by=list(class=counts$class, set=counts$Var2), FUN=sum)


p<-ggplot(counts_agg) +
  geom_bar(aes(x = class, y = x, fill = set), stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("ivory1", "deepskyblue2"), labels = c("non-expanded", "expanded")) +
  labs(x = "", y = "proportion", fill = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) 

svg("output_figures/Fig3B.svg")
print(p)
dev.off()

sharing2<-sharing[as.numeric(as.character(sharing$sharing_level)) >= 2,]

pX_0<-(786-sharing$sharing_level)/786 # this is equal to e^(-m)
m<--log(pX_0) # m is equal to -ln(pX_0) - where m is count in a repertoire of 10^5
f<-m*10 # average m in a repertoire of 10^6

sharing$f<-f
sharing[sharing$sharing_level<2,]$f<-min(sharing[sharing$sharing_level >= 2,]$f)/10 # for those that we cannot estimate, put a value a factor of 10 lower
sharing$f_permln<-sharing$f/10^6 # m per million

exp<-sharing[sharing$set=="exp",]
ctrl<-sharing[sharing$set=="ctrl",]

counts_exp <- data.frame(table(exp$f_permln))
counts_exp$log<-log10(as.numeric(as.character(counts_exp$Var1)))
counts_exp$logbin<-unlist(lapply(counts_exp$log, function(x){ifelse(x <= -8, "< 10^-8", 
                                                                    ifelse(x <= -7, "10^-8 - 10^-7", 
                                                                           ifelse(x <= -6, "10^-7 - 10^-6",
                                                                                  ifelse(x <= -5, "10^-6 - 10^-5", "> 10^-5"))))}))
counts_exp$logbin<-factor(counts_exp$logbin, 
                          levels = c("< 10^-8", "10^-8 - 10^-7", "10^-7 - 10^-6", "10^-6 - 10^-5", "> 10^-5"))
counts_exp_agg<-aggregate(counts_exp$Freq, by=list(logbin = counts_exp$logbin), FUN=sum)

p1<-ggplot(counts_exp_agg) +
  geom_bar(aes(x = logbin, y = x), stat = "identity", fill = "deepskyblue2", color = "black") +
  scale_x_discrete(labels = c(bquote('<'*10^-8),bquote(10^-8*' - '*10^-7), 
                              bquote(10^-7*' - '*10^-6), bquote(10^-6*' - '*10^-5), bquote('>'*10^-5))) +
  labs(x = "", y = "number of CDR3s", fill = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) 

svg("output_figures/Fig3C.svg")
print(p1)
dev.off()

counts_ctrl <- data.frame(table(ctrl$f_permln))
counts_ctrl$log<-log10(as.numeric(as.character(counts_ctrl$Var1)))
counts_ctrl$logbin<-unlist(lapply(counts_ctrl$log, function(x){ifelse(x <= -8, "< 10^-8", 
                                                                      ifelse(x <= -7, "10^-8 - 10^-7", 
                                                                             ifelse(x <= -6, "10^-7 - 10^-6",
                                                                                    ifelse(x <= -5, "10^-6 - 10^-5", "> 10^-5"))))}))
counts_ctrl$logbin<-factor(counts_ctrl$logbin, 
                           levels = c("< 10^-8", "10^-8 - 10^-7", "10^-7 - 10^-6", "10^-6 - 10^-5", "> 10^-5"))
counts_ctrl_agg<-aggregate(counts_ctrl$Freq, by=list(logbin = counts_ctrl$logbin), FUN=sum)

library(scales)
p2<-ggplot(counts_ctrl_agg) +
  geom_bar(aes(x = logbin, y = x), stat = "identity", fill = "ivory1", color = "black") +
  scale_x_discrete(labels = c(bquote('<'*10^-8),bquote(10^-8*' - '*10^-7), 
                              bquote(10^-7*' - '*10^-6), bquote(10^-6*' - '*10^-5), bquote('>'*10^-5))) +
  scale_y_continuous(labels = label_number()) +
  labs(x = "", y = "number of CDR3s", fill = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) 

svg("output_figures/Fig3D.svg")
print(p2)
dev.off()
