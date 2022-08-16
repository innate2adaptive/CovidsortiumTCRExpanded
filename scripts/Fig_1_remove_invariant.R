#Annotation of COVID-19 expanded TCRs

library(kernlab)
library(igraph)
library(stringdist)
library(tidyr)
library(dplyr)
library(data.table)
library(pheatmap)
#on Linux
drive <- "/media/benny/data/"
#on Windows
drive<-"C:/users/Benny Chain/"


#this is the path to the master files on my computer
input_data<-"Dropbox/Temp (1)/COVID-19/Data/"
input_ag_sp<-"Dropbox/Temp (1)/COVID-19/Antigen_specific/"
#dir<-dir(input_data)
#dir
#working folder for collecting output
#dir.create(paste0(drive,"Dropbox/R_temp/30_11_2021/"))
folder<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\data_for_paper\\"

file<-"exp_AB_wide1.RData"
file_name<-load(paste0(drive,input_data,file))

#remove MAIT
i_V<-which(exp_AB_wide1$v_call=="TRAV1-2")
i_J<-which(exp_AB_wide1$j_call=="TRAJ12" | exp_AB_wide1$j_call=="TRAJ20" |exp_AB_wide1$j_call=="TRAJ33")
inter<-intersect(i_V,i_J)
exp_AB_wide2<-exp_AB_wide1[-inter,]

#remove iKT
i_V<-which(exp_AB_wide1$v_call=="TRAV10")
i_J<-which(exp_AB_wide1$j_call=="TRAJ18")
exp_AB_wide3<-exp_AB_wide2[-(intersect(i_V,i_J)),]

save(exp_AB_wide3,file=  paste0(drive,input_data,"exp_AB_wide3.RData"))
