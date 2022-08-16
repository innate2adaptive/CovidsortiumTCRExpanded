#Annotation of COVID-19 expanded TCRs

library(kernlab)
library(igraph)
library(stringdist)
library(tidyr)
library(dplyr)
library(data.table)
library(pheatmap)


#this is the path to the master files on my computer
# input_data<-"Dropbox/Temp (1)/COVID-19/Data/"
# input_ag_sp<-"Dropbox/Temp (1)/COVID-19/Antigen_specific/"
# #dir<-dir(input_data)
#dir
#working folder for collecting output
#dir.create(paste0(drive,"Dropbox/R_temp/30_11_2021/"))

myURL<-"https://www.dropbox.com/s/zyz65xxlaq3u2x4/exp_AB_wide1.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

#remove MAIT
i_V<-which(exp_AB_wide1$v_call=="TRAV1-2")
i_J<-which(exp_AB_wide1$j_call=="TRAJ12" | exp_AB_wide1$j_call=="TRAJ20" |exp_AB_wide1$j_call=="TRAJ33")
inter<-intersect(i_V,i_J)
exp_AB_wide2<-exp_AB_wide1[-inter,]

#remove iKT
i_V<-which(exp_AB_wide1$v_call=="TRAV10")
i_J<-which(exp_AB_wide1$j_call=="TRAJ18")
exp_AB_wide3<-exp_AB_wide2[-(intersect(i_V,i_J)),]

save(exp_AB_wide3,file=  "data/output_data/exp_AB_wide3.RData")
