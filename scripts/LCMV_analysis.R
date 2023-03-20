library(ggplot2)
library(reshape)
library(stringr)

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

## Section A - total counts (Fig 5B)

counts_gp66<-data.frame(t(rbind(gp66naive_c,gp66eff_c)))
colnames(counts_gp66)<-c("naive", "effector")
counts_gp66$id<-rownames(counts_gp66)
counts_gp66[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_gp66$id, "[.]", 3)

g<-counts_gp66[(counts_gp66$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_gp66$naive))

counts_gp66_m<-melt(counts_gp66)

counts_gp92<-data.frame(t(rbind(gp92naive_c,gp92eff_c)))
colnames(counts_gp92)<-c("naive", "effector")
counts_gp92$id<-rownames(counts_gp92)
counts_gp92[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_gp92$id, "[.]", 3)

g<-counts_gp92[(counts_gp92$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_gp92$naive))

counts_gp92_m<-melt(counts_gp92)

counts_np205<-data.frame(t(rbind(np205naive_c,np205eff_c)))
colnames(counts_np205)<-c("naive", "effector")
counts_np205$id<-rownames(counts_np205)
counts_np205[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_np205$id, "[.]", 3)

g<-counts_np205[(counts_np205$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_np205$naive))

counts_np205_m<-melt(counts_np205)

counts_np396<-data.frame(t(rbind(np396naive_c,np396eff_c)))
colnames(counts_np396)<-c("naive", "effector")
counts_np396$id<-rownames(counts_np396)
counts_np396[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_np396$id, "[.]", 3)

g<-counts_np396[(counts_np396$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_np396$naive))

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