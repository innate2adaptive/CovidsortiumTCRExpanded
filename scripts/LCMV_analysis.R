library(ggplot2)
library(reshape)
library(stringr)

## Section A - total counts (Fig 5B)

gp66naive<-read.csv("data/LCMVdata/NaiveGP66_UMI.csv")
gp66eff<-read.csv("data/LCMVdata/E_GP66_UMI.csv")

gp92naive<-read.csv("data/LCMVdata/NaiveGP92_UMI.csv")
gp92eff<-read.csv("data/LCMVdata/E_GP92_UMI.csv")

np205naive<-read.csv("data/LCMVdata/NaiveNP205_UMI.csv")
np205eff<-read.csv("data/LCMVdata/E_NP205_UMI.csv")

np396naive<-read.csv("data/LCMVdata/NaiveNP396_UMI.csv")
np396eff<-read.csv("data/LCMVdata/E_NP396_UMI.csv")

gp66naive_c<-colSums(gp66naive[,2:11])
gp66eff_c<-colSums(gp66eff[,2:11])

gp92naive_c<-colSums(gp92naive[,2:11])
gp92eff_c<-colSums(gp92eff[,2:11])

np205naive_c<-colSums(np205naive[,2:11])
np205eff_c<-colSums(np205eff[,2:11])

np396naive_c<-colSums(np396naive[,2:11])
np396eff_c<-colSums(np396eff[,2:11])

counts_gp66<-data.frame(t(rbind(gp66naive_c,gp66eff_c)))
colnames(counts_gp66)<-c("naive", "effector")
counts_gp66$id<-rownames(counts_gp66)
counts_gp66[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_gp66$id, "[.]", 3)

counts_gp66_m<-melt(counts_gp66)

counts_gp92<-data.frame(t(rbind(gp92naive_c,gp92eff_c)))
colnames(counts_gp92)<-c("naive", "effector")
counts_gp92$id<-rownames(counts_gp92)
counts_gp92[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_gp92$id, "[.]", 3)

counts_gp92_m<-melt(counts_gp92)

counts_np205<-data.frame(t(rbind(np205naive_c,np205eff_c)))
colnames(counts_np205)<-c("naive", "effector")
counts_np205$id<-rownames(counts_np205)
counts_np205[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_np205$id, "[.]", 3)

counts_np205_m<-melt(counts_np205)

counts_np396<-data.frame(t(rbind(np396naive_c,np396eff_c)))
colnames(counts_np396)<-c("naive", "effector")
counts_np396$id<-rownames(counts_np396)
counts_np396[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_np396$id, "[.]", 3)

counts_np396_m<-melt(counts_np396)

counts_m<-rbind(counts_gp66_m, counts_gp92_m, counts_np205_m, counts_np396_m)


counts_m$LCMV<-factor(counts_m$LCMV, levels = c("PBS", "LCMV8", "LCMV40"))

p<-ggplot(counts_m) + 
  geom_point(aes(x = LCMV, y = value, color = variable), 
              position = position_dodge(width = 0.5), size = 5, alpha = 0.4) +
  labs(y = "TCR abundance", x = "", col = "") +
  scale_x_discrete(labels = c("PBS", "day 8", "day 40")) +
  scale_colour_manual(values = c(effector = "darkseagreen", naive = "darkmagenta")) +
  facet_wrap(vars(epitope), ncol=2) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=14),
        strip.text = element_text(size=16))

svg("output_figures/Fig5b_abundance_TCR_LCMV.svg")
print(p)
dev.off()

## Section B - compare direct and inference methods - This still does not work

# merge naive to measure sharing levels

gp66naive_m<-melt(gp66naive)
gp66naive_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(gp66naive_m$variable, "[.]", 3)

gp92naive_m<-melt(gp92naive)
gp92naive_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(gp92naive_m$variable, "[.]", 3)

np205naive_m<-melt(np205naive)
np205naive_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(np205naive_m$variable, "[.]", 3)

np396naive_m<-melt(np396naive)
np396naive_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(np396naive_m$variable, "[.]", 3)

naive<-rbind(gp66naive_m, gp92naive_m, np205naive_m, np396naive_m)
naive<-naive[naive$value > 0,]
naive[naive$LCMV == "PBS",]$mouse<- paste0(naive[naive$LCMV == "PBS",]$mouse, "a")
naive[naive$LCMV == "LCMV8",]$mouse<- paste0(naive[naive$LCMV == "LCMV8",]$mouse, "b")
naive[naive$LCMV == "LCMV40",]$mouse<- paste0(naive[naive$LCMV == "LCMV40",]$mouse, "c")


CDR3counts<-list()

for (cdr3 in unique(naive$X)){
  s<-naive[naive$X == cdr3,]
  print(s)
  mice<-unique(s$mouse)
  print(mice)
  print(length(mice))
  CDR3counts[[cdr3]]<-length(mice)
}

cdr3counts<-data.frame(t(data.frame(CDR3counts)))
colnames(cdr3counts)<-c("sharing_level")
cdr3counts$aminoAcid<-rownames(cdr3counts)

mice<-length(unique(naive$mouse))


# merge effectors

gp66eff_m<-melt(gp66eff)
gp66eff_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(gp66eff_m$variable, "[.]", 3)

gp92eff_m<-melt(gp92eff)
gp92eff_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(gp92eff_m$variable, "[.]", 3)

np205eff_m<-melt(np205eff)
np205eff_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(np205eff_m$variable, "[.]", 3)

np396eff_m<-melt(np396eff)
np396eff_m[c("epitope","LCMV", "mouse")]<-str_split_fixed(np396eff_m$variable, "[.]", 3)


effectors<-rbind(gp66eff_m, gp92eff_m, np205eff_m, np396eff_m)
effectors<-effectors[effectors$value>0,]

eff_sharing<-merge(effectors, cdr3counts[,c("aminoAcid", "sharing_level")], 
      by.x = "X", by.y = "aminoAcid", all.x=TRUE)
eff_sharing$sharing_level<-as.numeric(as.character(eff_sharing$sharing_level))

eff_sharing[is.na(eff_sharing$sharing_level),]<-0

pX_0<-(mice-eff_sharing$sharing_level)/mice # this is equal to e^(-m)
m<--log(pX_0) # m is equal to -ln(pX_0) - where m is count in a repertoire of 10^5
f<-m*10 # average m in a repertoire of 10^6

eff_sharing$f<-f
eff_sharing[eff_sharing$sharing_level<2,]$f<-min(eff_sharing[eff_sharing$sharing_level >= 2,]$f)/10 # for those that we cannot estimate, put a value a factor of 10 lower
eff_sharing$f_permln<-eff_sharing$f/10^6 # m per million

counts_exp <- data.frame(table(eff_sharing$f_permln))
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
  labs(x = "", y = "proportion", fill = "") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12))

print(p1)
