#Annotation of COVID-19 expanded TCRs

library(tidyr)
library(dplyr)
library(data.table)
library(pheatmap)

#on Linux
drive <- "/media/benny/data/"
#on Windows
drive<-"C:/users/Benny/"
#this is the path to the master files on my computer
input_data<-"Dropbox\\Temp (1)\\COVID-19\\data\\"
dir<-dir(paste0(drive,input_data))

#folder for plots
folder_plots<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Figs\\"
#folder for data output
folder_data<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Data_for_paper\\"

#dir.create(paste0(drive,"Dropbox/R_temp/30_11_2021/"))
#folder<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\data_for_paper\\"
#folder1<-"Dropbox/R_temp/30_11_2021/"
#folder_fig<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\figs\\"

#Section A - heatmap of samples and metadata

#load annotated TCR set
data<-read.table(paste0(drive, folder_data,"data.txt"),header=TRUE,sep="\t")
PCR<-read.table(paste0(drive, folder_data,"TCR_rep_HCW_18_12_2020_b.txt"),header=TRUE,sep="\t",stringsAsFactors = FALSE)
PCR_dat<-right_join(PCR,data,by =c("Sample.ID"="Sample_ID"))
rownames(data)<-data$Sample_ID
annotations<-data[,c(8,10,11)]
rownames(annotations)<-data$Sample_ID
mat_PCR<-PCR_dat[,2:7]
i_pos<-which(mat_PCR=="+++",arr.ind = TRUE)
mat_PCR[i_pos]<-"X"
i_NA<-which(mat_PCR!="X",arr.ind = TRUE)
mat_PCR[i_NA]<-""
i_NA<-which(is.na(mat_PCR),arr.ind = TRUE)
mat_PCR[i_NA]<-""
#dev.set(1)
pheatmap(data[,2:7],color = c("lightblue"),display_numbers=mat_PCR,fontsize_number=11,number_color="red",cluster_cols = FALSE,cluster_rows = FALSE,legend=FALSE, cellwidth = 10, cellheight = 10,  border_color="white",breaks = c(0,2,10),annotation_row = annotations)
imageSave(file=paste0(drive,folder_plots,"patient_heatmap.png"))
#save manually as svg Fig_1a


########################################################################################
#Fig 1c

file<-"exp_AB_wide3.RData"

file_name<-load(paste0(drive,folder_data,file))

i_ctrl<-which(exp_AB_wide3$control==T)

#calculate 
#extract controls
exp_AB_ctrl<-exp_AB_wide3[i_ctrl,]
a0<-unique(exp_AB_ctrl$total_0)
a1<-unique(exp_AB_ctrl$total_1)
a2<-unique(exp_AB_ctrl$total_2)
a3<-unique(exp_AB_ctrl$total_3)
a4<-unique(exp_AB_ctrl$total_4)

total_c<-sum(c(a0,a1,a2,a3,a4),na.rm=TRUE)
#calcualte false discovery rate. 
FDR<-dim(exp_AB_ctrl)[1]/total_c

i_chain<-which(exp_AB_wide3$chain=="beta")
#make vector contianing numbers correpnding to column names
col_names<-1:30
names(col_names)<-colnames(exp_AB_wide3)
col_names

#make file wiht unqiue ID and Control status
IDs<-distinct(data.frame(exp_AB_wide3$ID,exp_AB_wide3$control))
i_ctrl<-which(IDs[,2]=="FALSE")
#####################################################################
#Section B 
#plot a time course of total Abudnance at each time point. Individual HCW and box plot
#each chain separately
i_chain<-which(exp_AB_wide3$chain=="beta")
i_chain<-which(exp_AB_wide3$chain=="alpha")
PCR_total<-aggregate(exp_AB_wide3[i_chain,10:18], FUN = sum , by = list(exp_AB_wide3$ID[i_chain] ))
#alpha and beta together
PCR_total<-aggregate(exp_AB_wide3[,10:18], FUN = sum , by = list(exp_AB_wide3$ID ))

#replace 0 by -0.5
min<-sort(unique(unlist(PCR_total[,2:10])),decreasing = FALSE)[2]
#replace zeros by min/2
#i_zero<-which(PCR_total==0,arr.ind = TRUE)
#PCR_total[i_zero]<-min/2
#PCR_total$control<-IDs$exp_AB_wide1.control
#PCR_total_l<-pivot_longer(PCR_total,cols=2:10)
#PCR_total_l$log<-log2(PCR_total_l$value)

#analyse PCR+ first
PCR_total1<-PCR_total[i_ctrl,2:10]

#log
#plot(,0,xlim=c(0,10),ylim=c(0,log2(max(PCR_total1,na.rm=TRUE))),xaxt="none",las=1,xlab="", ylab="",yaxt="none")
#axis(1,at=c(1:9),label=c(-3:4,15),cex.axis=1.5)
#axis(2,at=c(log2(min/2),ceiling(log2(min/2)):15),labels=c(0,2^(ceiling(log2(min/2)):15)),cex.axis=1.5,las=1)
#i<-1
#for ( i in 1 : dim(PCR_total1)[1]){
#i_na<-which(!is.na(PCR_total1[i,]))
#x<-runif(i_na,i_na-(i_na/50),i_na+(i_na/50))
#y<-unlist(log2(PCR_total1[i,i_na]))
#points(x,y,pch=19,cex = 1.1,col="blue")
#}

#not logged looks nicer !
par(mar=c(5,4,4,2))
par(mar=c(3,2,2,2))
plot(-1,-1,xlim=c(0,10),ylim=c(0,(max(PCR_total1,na.rm=TRUE))),xaxt="none",las=1,xlab="Weeks from first PCR+", ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none",cex.lab=1.5)
#plot(-1,-1,xlim=c(0,10),ylim=c(0,80000),xaxt="none",las=1,xlab="Weeks to PCR+", ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none",cex.lab=1.5)

axis(1,at=c(1:9),label=c(-3:4,15),cex.axis=1.5)
axis(2,at=c(0,2,4,6,8,10)*1E4,labels=c(0,2,4,6,8,10),cex.axis=2,las=1)
i<-1
#c(log2(min/2),8:15)
for ( i in 1 : dim(PCR_total1)[1]){
  i_na<-which(!is.na(PCR_total1[i,]))
  x<-runif(i_na,i_na-(i_na/30),i_na+(i_na/30))
  y<-unlist(PCR_total1[i,i_na])
  points(x,y,pch=19,cex = 1.1,col="red")
}

boxplot((PCR_total1),add=TRUE,col=NA,axes=FALSE, at=1:9,border = "red",boxlwd=2,pars = list(boxwex = 0.4, staplewex = 1, outwex = 1),notch=F)

imageSave(file=paste0(drive,folder_plots,"timecourse-PCRpos.png"))


#add controls

i_ctrl1<-which(IDs[,2]=="TRUE")
PCR_total2<-PCR_total[i_ctrl1,5:9]

#not logged looks nicer !
plot(-1,-1,xlim=c(-1,5),ylim=c(0,(max(PCR_total1,na.rm=TRUE))),xaxt="none",las=1,xlab="Weeks from baseline", ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none",cex.lab=1.5)
#plot(-1,-1,xlim=c(0,10),ylim=c(0,80000),xaxt="none",las=1,xlab="Weeks to PCR+", ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none",cex.lab=1.5)

axis(1,at=c(0:5),label=c("BL","FUP1","FUP2","FUP3","FUP4",""),cex.axis=1.4)
axis(2,at=c(0,2,4,6,8,10)*1E4,labels=c(0,2,4,6,8,10),cex.axis=2.0,las=1)
i<-1
for ( i in 1 : dim(PCR_total2)[1]){
  i_na<-which(!is.na(PCR_total2[i,]))
  x<-(runif(i_na,i_na-(i_na/40),i_na+(i_na/40)))-1
  y<-unlist((PCR_total2[i,i_na]))
  points(x,y,pch=19,cex = 1.1,col="black")
}

boxplot((PCR_total2),add=TRUE,col=NA,axes=FALSE, boxlwd=2,at=c(0:4),border = "black",lwd=1,pars = list(boxwex = 0.4, staplewex = 1, outwex = 1),notch=F)
imageSave(file=paste0(drive,folder_plots,"timecourse-ctrl.png"))

#write.table(exp_AB_wide1,file=paste0(drive, input_data,"exp_AB_wide1.txt"),sep="\t")
###################################################
#correlation to symptoms
sum_zero<-PCR_total1$proportion_0
sum_all<-rowMeans(PCR_total1[,0:9],na.rm=TRUE)

data_merge<-merge(IDs,data,by.x= 1, by.y = 1)

#ethnicity
i_white<-which(data_merge$Ethnicity=="White")
i_black<-which(data_merge$Ethnicity=="Black")
i_asian<-which(data_merge$Ethnicity=="Asian")
i_bame<-union(i_black,i_asian)
ethnicity<-rep("blue",dim(data_merge)[1])
ethnicity[i_bame]<-"green"
#ethnicity
#y<-sum_zero
y<-sum_all
x<-data_merge$Symptoms
z<-lm(y~x)
plot(x,y,xlab="Symptoms",ylim = c(0,1E5), cex.lab=1.5,pch=19,xlim=c(0,15),xaxt="none",ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none")
axis(2,at=c(0,5,10,15,20,25,30)*1E4,labels=c(0,5,10,15,20,25,30),cex.axis=1.5,las=1)
axis(1, at = c(0,5,10,15),labels= c(0,5,10,15),cex.axis=1.5)
#?plot
z

abline(z,lwd=2)
cor.test(x,y,method = "spearman")
imageSave(file=paste0(drive, folder_fig,"symptoms.png"))

y_white<-y[i_white]
y_bame<-y[union(i_black,i_asian)]

wilcox.test(y_white,y_bame)

y1<-c(y_bame,y_white)
x1<-c(rnorm(length(y_bame),1,0.1),rnorm(length(y_white),2,0.1))
plot(x1,y1,pch=19,ylim=c(0,1E5),xlim=c(0,3),xaxt="none",xlab="",main=0.008,ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none",cex.lab=1.5)
axis(2,at=c(0,4,10,15,20,25,30)*1E4,labels=c(0,4,10,15,20,25,30),cex.axis=1.5,las=1)
axis(1,at=c(1,2),labels=c("BAME", "White"),cex.axis=1.5)
ethnicity <-rep(0,length(y))
ethnicity[i_white]<-1

boxplot(y~ethnicity,add=TRUE,col=NA,axes=FALSE,boxlwd=2,outline=FALSE)
summary(lm(y~ethnicity))
imageSave(file=paste0(drive, folder_fig,"Race.png"))

#sex
wilcox.test(data_merge$Symptoms[i_white],data_merge$Symptoms[union(i_black,i_asian)])

i_male<-which(data_merge$Sex=="Male")
i_female<-which(data_merge$Sex=="Female")

y2<-c(y[i_female],y[i_male])
x2<-c(rnorm(length(i_female),1,0.1),rnorm(length(i_male),2,0.1))
plot(x2,y2,pch=19,ylim=c(0,1E5),xlim=c(0,3),xaxt="none",xlab="",ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none",cex.lab=1.5)
axis(2,at=c(0,4,10,15,20,25,30)*1E4,labels=c(0,4,10,15,20,25,30),cex.axis=1.5,las=1)
axis(1,at=c(1,2),labels=c("Female", "Male"),cex.axis=1.5)

sex<-rep(0,length(y))
sex[i_male]<-1
wilcox.test(y[i_female],y[i_male])
boxplot(y~sex,add=TRUE,col=NA,axes=FALSE,boxlwd=2,outline=FALSE)
imageSave(file=paste0(drive, folder_fig,"gender.png"))
        
#####################
#Section D Serology time course
#NP antibody
data_merge<-merge(IDs,data,by.x= 1, by.y = 1)

#change week to week from PCR
NP<-read.table(file=paste0(drive,folder_data,"COVIDsortium_Roche_serology.txt"),check.names=FALSE,sep="\t",header=TRUE,stringsAsFactors = FALSE)
NP_merge1<-merge(data_merge,NP,by.x= 1, by.y = 1)
NP_merge<-NP_merge1[,14:32]
NP_long <-pivot_longer(NP_merge, cols = 4:19, names_to="week",values_to = "titre")
NP_long$weekfromPCR<-as.numeric(NP_long$week)-NP_long$Week_first
NP_long1<-NP_long[,-4]

#collect together all those samples >12 weeks 
#i_15<-which(NP_long1$weekfromPCR>12)
#NP_long1$weekfromPCR[i_15]<-15
NP_wide<-pivot_wider(NP_long1,names_from = weekfromPCR,values_from="titre")
NP_wide<-NP_wide[,c(1:3,22,4:21)]
NP_titre<-NP_wide[,4:22]
#average weeks 12:15
NP_15<-rowMeans(NP_wide[,18:22],na.rm=TRUE)
NP_titre_s<-NP_titre[,1:9]
NP_titre_s[,10]<-NP_15
#plot
plot(-1,-1,xlim=c(0,10),ylim=c(0,(max(NP_titre_s,na.rm=TRUE))),xaxt="none",las=1,xlab="Weeks to PCR+", ylab="NP Antibody titre",yaxt="none",cex.lab=1.5)
axis(1,at=c(1:10),label=c(-3:5,15),cex.axis=1.5)
axis(2,at=c(0,50,100,150),labels=c(0,50,100,150),cex.axis=1.4,las=1)
i<-1
#c(log2(min/2),8:15)
for ( i in 1 : dim(NP_titre)[1]){
  i_na<-which(!is.na(NP_titre_s[i,]))
  x<-runif(i_na,i_na-(i_na/30),i_na+(i_na/30))
  y<-unlist(NP_titre_s[i,i_na])
  points(x,y,pch=19,cex = 1.1,col="blue",type = "p")
}

boxplot((NP_titre_s),add=TRUE,col=NA,axes=FALSE, at=1:10,border = "blue",boxlwd=2,pars = list(boxwex = 0.4, staplewex = 1, outwex = 1),notch=F)
imageSave(file=paste0(drive, folder_plots,"NP_timecourse.png"))
#dir(paste0(drive,folder1))
early_NP<-rowMeans(NP_titre[,1:8],na.rm=TRUE)
total_NP<-rowMeans(NP_titre[,1:19],na.rm=TRUE)

plot(early_NP,total_NP)

plot(early_NP,sum_all,xlab="NP Antibody",ylab="Total TCR abundance",pch=19,col="blue",cex.axis=1.5,cex.lab=1.5,main=0.019)
abline(lm(sum_all~early_NP),lwd=2)
imageSave(file=paste0(drive, folder_fig,"NP_corr.png"))

cor.test(total_NP,sum_all,method="spearman")

cor.test(early_NP,sum_all,method="spearman")


#############################################################

#SP antibody
#change week to week from PCR
NP<-read.table(file=paste0(drive,folder_data,"Serology_clinical_09_2020_SP.csv"),sep=",",na.strings="",check.names=FALSE,header=TRUE,stringsAsFactors = FALSE)
NP_merge1<-merge(data_merge,NP[1:500,],by.x= 1, by.y = 1)
NP_merge<-NP_merge1[,14:31]

NP_long <-pivot_longer(NP_merge, cols = 4:18, names_to="week",values_to = "titre",names_transform = list("titre" = as.numeric))
NP_long$titre<-as.numeric(NP_long$titre)
NP_long$titre[1:10]
NP_long$weekfromPCR<-as.numeric(NP_long$week)-NP_long$Week_first
NP_long1<-NP_long[,-4]
#is.numeric(unique(NP_long1$weekfromPCR))
#collect together all those samples >12 weeks 
#i_15<-which(NP_long1$weekfromPCR>12)
#NP_long1$weekfromPCR[i_15]<-15
NP_wide2<-pivot_wider(NP_long1,names_from = weekfromPCR,values_from="titre",names_repair = "minimal")
NP_wide1<-NP_wide2[,-3]
#NP_wide1[4,]
#?pivot_wider
NP_wide<-NP_wide1[,c(1:2,20,3:19)]
NP_titre<-NP_wide[,3:20]
NP_titre[3,]
#NP_titre[1,]
#average weeks 12:15
NP_15<-rowMeans(NP_wide[,18:20],na.rm=TRUE)

#NP_titre_s<-NP_titre[,1:17]
#NP_titre_s[,18]<-NP_15

NP_titre_s<-NP_titre[,1:8]
NP_titre_s[,9]<-NP_15


i<-1
for ( i in 1:10){NP_titre_s[,i]<-as.numeric(unlist(NP_titre_s[,i]))}

#plot
#plot(-1,-1,xlim=c(0,10),ylim=c(0,max(NP_titre_s,na.rm = TRUE)),xaxt="none",las=1,xlab="Weeks to PCR+", ylab="SP Antibody titre",yaxt="none",cex.lab=1.5)
plot(-1,-1,xlim=c(0,10),ylim=c(0,12),xaxt="none",las=1,xlab="Weeks to PCR+", ylab="SP Antibody titre",yaxt="none",cex.lab=1.5)

axis(1,at=c(1:10),label=c(-3:5,15),cex.axis=1.3)
axis(2,at=c(0,5,10,15,20),labels=c(0,5,10,15,20),cex.axis=1.4,las=1)
i<-1
#c(log2(min/2),8:15)
for ( i in 1 : dim(NP_titre_s)[1]){
  i_na<-which(!is.na(NP_titre_s[i,]))
  x<-runif(i_na,i_na-(i_na/50),i_na+(i_na/50))+1
  y<-as.numeric(unlist(NP_titre_s[i,i_na]))
  points(x,y,pch=19,cex = 1.1,col="blue",type = "p")
  i
}

a<-boxplot(NP_titre_s,add=TRUE,col=NA,axes=FALSE, at=2:10,border = "blue",boxlwd=2,pars = list(boxwex = 0.4, staplewex = 1, outwex = 1),notch=F)
imageSave(file=paste0(drive, folder_fig,"SP_time.png"))
#dir(paste0(drive,folder1))

early_SP<-rowMeans(NP_titre[,1:8],na.rm=TRUE)
total_SP<-rowMeans(NP_titre[,1:18],na.rm=TRUE)

plot(early_SP,total_SP)

plot(early_SP,sum_all,xlab="SP Antibody",ylab="Total TCR abundance",pch=19,col="blue",cex.axis=1.5,cex.lab=1.5)
abline(lm(sum_all~early_SP),lwd=2)
imageSave(file=paste0(drive, folder_fig,"SP_corr.png"))

cor.test(total_SP,sum_all,method="spearman")

cor.test(early_SP,sum_all,method="spearman")

#summary file
summary_table1<-data.frame(cbind(data_merge,sum_all,sum_zero,total_NP,early_NP,total_SP,early_SP))
summary_table<-summary_table1[,c(1,9,11:21)]

write.table(summary_table,file=paste0(drive,folder,"summary_table.txt"),sep="\t",row.names = FALSE)

