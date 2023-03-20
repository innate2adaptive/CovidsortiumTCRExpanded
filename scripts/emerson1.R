# load and compute emerson sharing

library(ggplot2)
library(reshape)
library(dplyr)
library(data.table)

options(timeout = max(10000, getOption("timeout"))) # increasing timeout might be necessary to load the files

myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

# the data is in long format, so first I keep only one entry for each tcr in each patient (here I have multiple entries for multiple timepoints)
all_B_long1<-all_B_long[!(duplicated(all_B_long[, c("decombinator_id", "ID")])),]
rm(all_B_long)
all_B_long1<-all_B_long1[!duplicated(all_B_long1[,c("v_call", "junction_aa")]),]

# load emerson
# sharing_ems<-readRDS("data/downloaded_data/Emerson_sharing_aa.rds")
sharing_ems1<-fread("data/output_data/emerson_publicity_counts.csv.gz")
colnames(sharing_ems1)<-c("V1", "amino_acid", "v_gene", "sharing_level")
# load expanded
load("data/output_data/exp_AB_wide4.RData")

exp_b<-exp_AB_wide3[exp_AB_wide3$chain == "beta",] # emerson set only has beta
exp_b<-exp_b[!duplicated(exp_b[,c("v_call", "junction_aa")]),]
exp_b_early<-exp_b[exp_b$max_timepoint_class == 'early',]
dim(exp_b_early)
exp_b_late<-exp_b[exp_b$max_timepoint_class == 'late',]
dim(exp_b_late)
exp_b_early<-exp_b_early[!duplicated(exp_b_early[,c("v_call", "junction_aa")]),]
exp_b_late<-exp_b_late[!duplicated(exp_b_late[,c("v_call", "junction_aa")]),]

rm(exp_AB_wide3)

# merge with sharing - adds a column which says how many emerson individuals the cdr3 is found in
sharing_exp<-merge(exp_b, sharing_ems1, 
                   by.x = c("v_call", "junction_aa"), 
                   by.y = c("v_gene", "amino_acid"), 
                   all.x = TRUE) # add sharing info
sharing_exp[is.na(sharing_exp$sharing_level),]$sharing_level<-0

sharing_exp_early<-merge(exp_b_early, sharing_ems1, 
                         by.x = c("v_call", "junction_aa"), 
                         by.y = c("v_gene", "amino_acid"), 
                         all.x = TRUE) # add sharing info
sharing_exp_early[is.na(sharing_exp_early$sharing_level),]$sharing_level<-0

sharing_exp_late<-merge(exp_b_late, sharing_ems1, 
                        by.x = c("v_call", "junction_aa"), 
                        by.y = c("v_gene", "amino_acid"), 
                        all.x = TRUE) # add sharing info
sharing_exp_late[is.na(sharing_exp_late$sharing_level),]$sharing_level<-0

dim(sharing_exp[sharing_exp$sharing_level >= 2,]) # 928 as in paper
dim(sharing_exp_early[sharing_exp_early$sharing_level >= 2,]) # 135
dim(sharing_exp_late[sharing_exp_late$sharing_level >= 2,]) # 71

rm(exp_b)
rm(exp_b_early)
rm(exp_b_late)

all_B_long1$ID<-as.character(all_B_long1$ID)
sharing_exp$ID<-as.character(sharing_exp$ID)
nonexp_B_long<-anti_join(all_B_long1, sharing_exp, by=c("decombinator_id", "ID")) # remove expanded from the list
rm(all_B_long1)

# merge with sharing - adds a column which says how many emerson individuals the cdr3 is found in
sharing_ctrl_B<-merge(nonexp_B_long, sharing_ems1, 
                      by.x = c("v_call", "junction_aa"), 
                      by.y = c("v_gene", "amino_acid"), 
                      all.x = TRUE) # add sharing info
sharing_ctrl_B[is.na(sharing_ctrl_B$sharing_level),]$sharing_level<-0
dim(sharing_ctrl_B[sharing_ctrl_B$sharing_level >= 2,]) # 1136870

sharing_ctrl_B$set<-"ctrl"
sharing_exp$set<-"exp"
sharing_exp_early$set<-"early"
sharing_exp_late$set<-"late"

sharing<-rbind(sharing_ctrl_B[c("junction_aa", "v_call", "sharing_level", "ID", "decombinator_id", "set")], 
               sharing_exp[c("junction_aa", "v_call", "sharing_level", "ID", "decombinator_id", "set")],
               sharing_exp_early[c("junction_aa", "v_call", "sharing_level", "ID", "decombinator_id", "set")],
               sharing_exp_late[c("junction_aa", "v_call", "sharing_level", "ID", "decombinator_id", "set")]
               )

counts<-data.frame(table(sharing$sharing_level, sharing$set))
counts$class <- unlist(lapply(as.numeric(as.character(counts$Var1)), function(x) {ifelse(x <= 3, "0-3",
                                                               ifelse((x>3) & (x <= 10), "4-10", 
                                                                      ifelse((x>10) & (x <= 30), "11-30",
                                                                             ifelse((x>30) & (x <= 100), "31-100", 
                                                                                    ifelse((x>100) & (x <= 300), "101-300", "301+")))))}))

sum_ctrl <- sum(as.numeric(as.character(counts[counts$Var2 == "ctrl",]$Freq)))
sum_exp <- sum(as.numeric(as.character(counts[counts$Var2 == "exp",]$Freq)))
sum_exp_early <- sum(as.numeric(as.character(counts[counts$Var2 == "early",]$Freq)))
sum_exp_late <- sum(as.numeric(as.character(counts[counts$Var2 == "late",]$Freq)))

counts$prop<-0

counts[counts$Var2 == "ctrl",]$prop <- counts[counts$Var2 == "ctrl",]$Freq/sum_ctrl
counts[counts$Var2 == "exp",]$prop <- counts[counts$Var2 == "exp",]$Freq/sum_exp
counts[counts$Var2 == "early",]$prop <- counts[counts$Var2 == "early",]$Freq/sum_exp_early
counts[counts$Var2 == "late",]$prop <- counts[counts$Var2 == "late",]$Freq/sum_exp_late

counts$class <- factor(counts$class, levels = c("0-3", "4-10", "11-30", "31-100", "101-300", ">301"))
counts$Var2 <- factor(counts$Var2, levels = c("ctrl", "exp", "early", "late"))


counts_agg<-aggregate(counts$prop, by=list(class=counts$class, set=counts$Var2), FUN=sum)


p<-ggplot(counts_agg) +
  geom_bar(aes(x = class, y = x, fill = set), stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("ivory1", "deepskyblue2", "orange", "cyan2"))+#, labels = c("non-expanded", "expanded")) +
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

write.csv(sharing, "data/output_data/Emerson_sharing_levels_calculated.csv")

exp<-sharing[sharing$set=="exp",]
print(mean(exp$f))
early<-sharing[sharing$set=="early",]
late<-sharing[sharing$set=="late",]
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

counts_early <- data.frame(table(early$f_permln))
counts_early$log<-log10(as.numeric(as.character(counts_early$Var1)))
counts_early$logbin<-unlist(lapply(counts_early$log, function(x){ifelse(x <= -8, "< 10^-8", 
                                                                    ifelse(x <= -7, "10^-8 - 10^-7", 
                                                                           ifelse(x <= -6, "10^-7 - 10^-6",
                                                                                  ifelse(x <= -5, "10^-6 - 10^-5", "> 10^-5"))))}))
counts_early$logbin<-factor(counts_early$logbin, 
                          levels = c("< 10^-8", "10^-8 - 10^-7", "10^-7 - 10^-6", "10^-6 - 10^-5", "> 10^-5"))
counts_early_agg<-aggregate(counts_early$Freq, by=list(logbin = counts_early$logbin), FUN=sum)

p1<-ggplot(counts_early_agg) +
  geom_bar(aes(x = logbin, y = x), stat = "identity", fill = "darkorange", color = "black") +
  scale_x_discrete(labels = c(bquote('<'*10^-8),bquote(10^-8*' - '*10^-7), 
                              bquote(10^-7*' - '*10^-6), bquote(10^-6*' - '*10^-5), bquote('>'*10^-5))) +
  labs(x = "", y = "number of CDR3s", fill = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) 

svg("output_figures/Fig3C_early.svg")
print(p1)
dev.off()

counts_late <- data.frame(table(late$f_permln))
counts_late$log<-log10(as.numeric(as.character(counts_late$Var1)))
counts_late$logbin<-unlist(lapply(counts_late$log, function(x){ifelse(x <= -8, "< 10^-8", 
                                                                    ifelse(x <= -7, "10^-8 - 10^-7", 
                                                                           ifelse(x <= -6, "10^-7 - 10^-6",
                                                                                  ifelse(x <= -5, "10^-6 - 10^-5", "> 10^-5"))))}))
counts_late$logbin<-factor(counts_late$logbin, 
                          levels = c("< 10^-8", "10^-8 - 10^-7", "10^-7 - 10^-6", "10^-6 - 10^-5", "> 10^-5"))
counts_late_agg<-aggregate(counts_late$Freq, by=list(logbin = counts_late$logbin), FUN=sum)

p1<-ggplot(counts_late_agg) +
  geom_bar(aes(x = logbin, y = x), stat = "identity", fill = "cyan2", color = "black") +
  scale_x_discrete(labels = c(bquote('<'*10^-8),bquote(10^-8*' - '*10^-7), 
                              bquote(10^-7*' - '*10^-6), bquote(10^-6*' - '*10^-5), bquote('>'*10^-5))) +
  labs(x = "", y = "number of CDR3s", fill = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) 

svg("output_figures/Fig3C_late.svg")
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
