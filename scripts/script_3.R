#General analysis of number of expanded TCRs and specifically 017
library(gplots)
library(RColorBrewer)
library(kernlab)
library(igraph)
library(stringdist)
library(ggplot2)

#folder for plots
folder_plots<-"output_figures/" # don't save to dropbox, so if people use this will save locally
#folder for data output
folder_data<-"data/output_data/"  # don't save to dropbox, so if people use this will save locally

file<-"exp_AB_wide3.RData"

file_name<-load(paste0(folder_data,file))

exp_AB<-get(file_name)
counts_all<-as.matrix(exp_AB[,10:18])
IDs<-unique(exp_AB$ID)
max<-max(counts_all,na.rm=TRUE)


results<-data.frame()

i<-1
for ( i in 1:length(IDs)){
#ID<-"0017"
#ID<-"0042"
#ID<-"0084"
#ID<-"0101"
#ID<-"0036"
#ID<-"0123"
ID<-IDs[i]
print(ID)
i_ID<-which(exp_AB$ID==ID)

counts<-counts_all[i_ID,]
chain<-exp_AB$chain[i_ID]
control<-exp_AB$control[i_ID]
#i_zero<-which(counts==0,arr.ind = TRUE)
min<-unique(sort(counts))[2]/2
counts_l<-counts
#counts_l<-log2(counts)
#counts_l[i_zero]<-min
max<-max(counts_l,na.rm=TRUE)
#counts<-as.matrix(exp_AB[i_17,13:17])
####################################################
#################################################################

# work with counts df so that I can then plot with ggplot2 at the end

counts_df <- data.frame(counts_l)
counts_df["tcrname"]<-rownames(counts_df)
counts_df["chain"]<-chain
counts_df_m<-reshape2::melt(counts_df)
counts_df_m["ID"]<-ID
counts_df_m["control"]<-unique(control)

results<-rbind(results, counts_df_m)

}


# Fig S8A
p<-ggplot(results[(results$control == FALSE) & (!is.na(results$value)),]) +
  geom_line(aes(x = variable, y = value, group = tcrname, col = chain)) +
  geom_point(aes(x = variable, y = value, col = chain)) +
  scale_color_manual(values = c(alpha = "brown2", beta="navyblue")) +
  scale_x_discrete(breaks = c("proportion_.3", "proportion_.2", "proportion_.1", 
                              "proportion_0", "proportion_1", "proportion_2", "proportion_3", "proportion_4", "proportion_14"),
                   labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "12-14")) +
  scale_y_continuous(limits = c(0,6000)) +
  labs(x = "weeks", y = "TCR/million") +
  facet_wrap(vars(ID), ncol = 7) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=16),
        title=element_text(size=14))

# Fig S8B
p_c<-ggplot(results[(results$control == TRUE) & (!is.na(results$value)),]) +
  geom_line(aes(x = variable, y = value, group = tcrname, col = chain)) +
  geom_point(aes(x = variable, y = value, col = chain)) +
  scale_color_manual(values = c(alpha = "brown2", beta="navyblue")) +
  scale_x_discrete(breaks = c("proportion_.3", "proportion_.2", "proportion_.1", 
                              "proportion_0", "proportion_1", "proportion_2", "proportion_3", "proportion_4", "proportion_14"),
                   labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "12-14")) +
  scale_y_continuous(limits = c(0,6000)) +
  labs(x = "weeks", y = "TCR/million") +
  facet_wrap(vars(ID)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=16),
        title=element_text(size=14))

svg(paste0(folder_plots, "FigS8A_dynamics_PCR+.svg"), width = 20, height = 15)
print(p)
dev.off()

svg(paste0(folder_plots, "FigS8B_dynamics_PCR-.svg"))
print(p_c)
dev.off()

# Fig 1E

p1<-ggplot(results[(results$ID == "0101"),]) +
  scale_x_discrete(breaks = c("proportion_.3", "proportion_.2", "proportion_.1", 
                              "proportion_0", "proportion_1", "proportion_2", "proportion_3", "proportion_4", "proportion_14"),
                   labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "12-14")) +
  scale_y_continuous(limits = c(0,600)) +
  geom_point(aes(x = variable, y = value, col = chain)) +
  geom_line(aes(x = variable, y = value, group = tcrname, col = chain), 
            data = results[(results$ID == "0101") & (!is.na(results$value)),]) +
  scale_color_manual(values = c(alpha = "brown2", beta="navyblue")) +
  labs(x = "weeks", y = "TCR/million") +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=16),
        title=element_text(size=14))

svg(paste0(folder_plots, "Fig1E_dynamics_HCW101.svg"))
print(p1)
dev.off()

# Fig S7

p1a<-ggplot(results[(results$ID == "0101") & (results$variable %in% c("proportion_.2", "proportion_14")), ]) +
  scale_x_discrete(breaks = c("proportion_.3", "proportion_14"),
                   labels = c("Baseline", "Week 12-14")) +
  # scale_y_continuous(limits = c(0,500)) +
  # geom_point(aes(x = variable, y = value, col = chain)) +
  geom_line(aes(x = variable, y = value, group = tcrname, col = chain), 
            data = results[(results$ID == "0101") & (!is.na(results$value)) & (results$variable %in% c("proportion_.3", "proportion_14")),]) +
  scale_color_manual(values = c(alpha = "brown2", beta="navyblue")) +
  labs(x = "weeks", y = "TCR/million") +
  theme_classic() + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=16),
        title=element_text(size=14))

svg(paste0(folder_plots, "FigS7_dynamics_HCW101_week-2_14.svg"))
print(p1a)
dev.off()

X<-results[(results$ID == "0101") & (results$variable %in% c("proportion_.3", "proportion_14")),]$value
# X[is.na(X)]<-0
kruskal.test(X,
             results[(results$ID == "0101") & (results$variable %in% c("proportion_.2", "proportion_14")), ]$variable)

# Fig 1F

p2<-ggplot(results[(results$ID == "0017"),]) +
  scale_x_discrete(breaks = c("proportion_.3", "proportion_.2", "proportion_.1", 
                              "proportion_0", "proportion_1", "proportion_2", "proportion_3", "proportion_4", "proportion_14"),
                   labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "12-14")) +
  scale_y_continuous(limits = c(0,600)) +
  geom_point(aes(x = variable, y = value, col = chain)) +
  geom_line(aes(x = variable, y = value, group = tcrname, col = chain), 
            data = results[(results$ID == "0017") & (!is.na(results$value)),]) +
  scale_color_manual(values = c(alpha = "brown2", beta="navyblue")) +
  labs(x = "weeks", y = "TCR/million") +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=16),
        title=element_text(size=14))

svg(paste0(folder_plots, "Fig1F_dynamics_HCW17.svg"))
print(p2)
dev.off()
