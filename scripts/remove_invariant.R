#Annotation of COVID-19 expanded TCRs

library(kernlab)
library(igraph)
library(stringdist)
library(tidyr)
library(dplyr)
library(data.table)
library(pheatmap)

load("data/output_data/exp_AB_wide1.RData")

#remove MAIT
i_V<-which(exp_AB_wide1$v_call=="TRAV1-2")
i_J<-which(exp_AB_wide1$j_call=="TRAJ12" | exp_AB_wide1$j_call=="TRAJ20" |exp_AB_wide1$j_call=="TRAJ33")
inter<-intersect(i_V,i_J)
exp_AB_wide2<-exp_AB_wide1[-inter,]

#remove iKT
i_V<-which(exp_AB_wide1$v_call=="TRAV10")
i_J<-which(exp_AB_wide1$j_call=="TRAJ18")
exp_AB_wide3<-exp_AB_wide2[-(intersect(i_V,i_J)),]


exp_AB_wide3$ID<-as.character(exp_AB_wide3$ID)
save(exp_AB_wide3,file=  "data/output_data/exp_AB_wide3.RData")
