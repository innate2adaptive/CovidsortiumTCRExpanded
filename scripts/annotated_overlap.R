## Fig 2A - show numbers of annotated CDR3 in each published set

load("data/VDJdb_single_specificity.RData")

vdj.cov<-VDJdb_all2[!is.na(VDJdb_all2$`SARS-CoV-2`),]
length(unique(vdj.cov$CDR3))

a.vdjdb<-vdj.cov[vdj.cov$Gene == "TRA",]
b.vdjdb<-vdj.cov[vdj.cov$Gene == "TRB",]
length(unique(a.vdjdb$CDR3))
length(unique(b.vdjdb$CDR3))

tao<-read.csv("data/Tao_CDR3.txt", sep = "\t")
length(unique(tao$CDR3))
a.tao<-tao[tao$gene == "alpha",]
b.tao<-tao[tao$gene == "beta",]
length(unique(a.tao$CDR3))
length(unique(b.tao$CDR3))

francis<-read.csv("data/Francis_Science_CDR3.txt", sep = "\t")
length(unique(francis$CDR3))
a.francis<-francis[francis$chain == "alpha",]
b.francis<-francis[francis$chain == "beta",]
length(unique(a.francis$CDR3))
length(unique(b.francis$CDR3))


length(c(unique(francis$CDR3), unique(tao$CDR3), unique(vdj.cov$CDR3)))

## Fig 2B - Venn diagram of overlap

load("data/output_data/exp_AB_wide3.RData")
a.exp.nc<-exp_AB_wide3[(exp_AB_wide3$chain == "alpha") & (exp_AB_wide3$control==FALSE),]
b.exp.nc<-exp_AB_wide3[(exp_AB_wide3$chain == "beta") & (exp_AB_wide3$control==FALSE),]
a.exp.c<-exp_AB_wide3[(exp_AB_wide3$chain == "alpha") & (exp_AB_wide3$control==TRUE),]
b.exp.c<-exp_AB_wide3[(exp_AB_wide3$chain == "beta") & (exp_AB_wide3$control==TRUE),]

library(ggvenn)

svg("output_figures/Fig2B_VennDiagram.svg")
ggvenn(
  list(Annotated = c(francis$CDR3, tao$CDR3, vdj.cov$CDR3), 
       Expanded = c(a.exp.nc$junction_aa, b.exp.nc$junction_aa)), 
  show_percentage = FALSE,
  fill_color = c("turquoise2", "palegreen1"),
  stroke_color = "lightgrey",
  stroke_size = 1, set_name_size = 8, 
  text_size = 8,
)
dev.off()

# number in PCR+ - non-unique
x<-sum(c(a.exp.nc$junction_aa, b.exp.nc$junction_aa) %in% c(francis$CDR3, tao$CDR3, vdj.cov$CDR3))
print(x)
x/(length(unique(c(francis$CDR3, tao$CDR3, vdj.cov$CDR3))))*100 # overlap over total number of annotated

# number in PCR- controls - non-unique
x0<-sum(unique(c(a.exp.c$junction_aa, b.exp.c$junction_aa)) %in% c(francis$CDR3, tao$CDR3, vdj.cov$CDR3))
print(x0)

# load controls and compare the numbers of annotated
# check this runs with new controls

load("data/output_data/control_alpha_1.RData")
load("data/output_data/control_beta_1.RData")

overlap <-list()

for (i in 1:10){
  a<-sample(control_a[[i]], length(a.exp.nc$junction_aa))
  b<-sample(control_b[[i]], length(b.exp.nc$junction_aa))
  
  x1<-sum(unlist(c(a, b)) %in% c(francis$CDR3, tao$CDR3, vdj.cov$CDR3))
  overlap<-c(overlap, x1)
}

mean(unlist(overlap))
mean(unlist(overlap))/(length(unique(c(francis$CDR3, tao$CDR3, vdj.cov$CDR3))))*100 # overlap over total number of annotated

# bar plots - look at unique intersect (Benny was using non-unique)
# need to sum non-unique and divide by correct chain to get same results as Bennt
all.cov<-c(francis$CDR3, tao$CDR3, vdj.cov$CDR3)
vdj.hiv<-VDJdb_all2[!is.na(VDJdb_all2$`HIV-1`),]
vdj.ebv<-VDJdb_all2[!is.na(VDJdb_all2$EBV),]
vdj.cmv<-VDJdb_all2[!is.na(VDJdb_all2$CMV),]

x.cov.a<-sum(a.exp.nc$junction_aa %in% all.cov)/length(all.cov)
print(x.cov.a)
x.cov.b<-sum(b.exp.nc$junction_aa %in% all.cov)/length(all.cov)
print(x.cov.b)

a.vdj.cmv<-vdj.cmv[vdj.cmv$Gene == "TRA",]
x.cmv.a<-length(intersect(unique(a.exp.nc$junction_aa), a.vdj.cmv$CDR3))/length(a.vdj.cmv$CDR3)
print(x.cmv.a)
x.cmv.b<-length(intersect(unique(b.exp.nc$junction_aa), vdj.cmv$CDR3))/length(vdj.cmv$CDR3)
print(x.cmv.b)

x.ebv.a<-length(intersect(unique(a.exp.nc$junction_aa), vdj.ebv$CDR3))/length(vdj.ebv$CDR3)
print(x.ebv.a)
x.ebv.b<-length(intersect(unique(b.exp.nc$junction_aa), vdj.ebv$CDR3))/length(vdj.ebv$CDR3)
print(x.ebv.b)

x.hiv.a<-length(intersect(unique(a.exp.nc$junction_aa), vdj.hiv$CDR3))/length(vdj.hiv$CDR3)
print(x.hiv.a)
x.hiv.b<-length(intersect(unique(b.exp.nc$junction_aa), vdj.hiv$CDR3))/length(vdj.hiv$CDR3)
print(x.hiv.b)

results<-data.frame()
results["SARS-CoV-2", "alpha"]<-x.cov.a
results["SARS-CoV-2", "beta"]<-x.cov.b
results["CMV", "alpha"]<-x.cmv.a
results["CMV", "beta"]<-x.cmv.b
results["EBV", "alpha"]<-x.ebv.a
results["EBV", "beta"]<-x.ebv.b
results["HIV", "alpha"]<-x.hiv.a
results["HIV", "beta"]<-x.hiv.b
results["virus"]<-rownames(results)

library(reshape)
results_m<-reshape::melt(results)

results_m$virus<-factor(results_m$virus, levels = c("SARS-CoV-2", "CMV", "EBV", "HIV"))
results_m$variable<-factor(results_m$variable, levels = c("beta", "alpha"))

library(ggplot2)

ggplot(results_m) +
  geom_bar(aes(x = virus, y = value, fill = variable), stat = "identity")
