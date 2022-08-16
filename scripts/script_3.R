#General analysis of number of expanded TCRs and specifically 017
library(gplots)
library(RColorBrewer)
library(kernlab)
library(igraph)
library(stringdist)

#on Linux
drive <- "/media/benny/data/"
#on Windows
drive<-"C:/users/Benny/"

#folder for plots
folder_plots<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Figs\\"
#folder for data output
folder_data<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Data_for_paper\\"

file<-"exp_AB_wide3.RData"

file_name<-load(paste0(drive,folder_data,file))

exp_AB<-get(file_name)
counts_all<-as.matrix(exp_AB[,10:18])
IDs<-unique(exp_AB$ID)
max<-max(counts_all,na.rm=TRUE)
i<-1
for ( i in 1:length(IDs)){
#ID<-"0017"
#ID<-"0042"
#ID<-"0084"
#ID<-"0195"
#ID<-"0036"
#ID<-"0123"
ID<-IDs[i]
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
#plotting individual traces
#exp_ID<-exp_AB_wide_pos1[i_17,]
#plot(-1,-1, xlim = c(1,9), ylim = c(log2(min),max(counts_l,na.rm=TRUE)),main=paste("Cluster ",c))
if (control[1]=="TRUE") {
plot(-2,-2,xlim=c(0,6),ylim=c(0,6000),xaxt="n",yaxt="n",xlab= "weeks", ylab="TCR abundance (per 10^3)",main=paste(ID),cex.lab=1.5)
axis(1,at = c(1,2,3,4,5),labels = c("BL", "FUP1","FUP2","FUP3","FUP4"),cex.axis=1.5)
axis(2,las=2,at = c(0,2000,4000,6000),labels=c(0,2,4,6),cex.axis=1.5)
}

#plot(-1,-1, xlim = c(1,9), ylim = c(log2(min),max(counts_l,na.rm=TRUE)),main=paste("Cluster ",c))
if (control[1]=="FALSE") {
  plot(-2,-2,xlim=c(0,10),ylim=c(0,6000),xaxt="n",yaxt="n",xlab= "weeks from PCR+",ylab="TCR abundance (per 10^3)",main=paste(ID),cex.lab=1.5)
  axis(1,at = c(1,2,3,4,5,6,7,8,9),labels = c(-3,-2,-1,0,1,2,3,4,14),cex.axis=1.5)
  axis(2,las=2,at = c(0,2000,4000,6000),labels=c(0,2,4,6),cex.axis=2.0)
  }

j<-1
for (j in 1: length(i_ID)){
    x1<-which(!is.na(counts_l[j,]))
    y<-rnorm(length(counts_l[j,x1]),counts_l[j,x1],counts_l[j,x1]/20)
    if (control[1]=="TRUE"){  x<-rnorm(length(x1),x1,x1/50) -3}
    if (control[1]=="FALSE"){  x<-rnorm(length(x1),x1,x1/50)}
    if (chain[j]=="alpha"){  points(x,y,"b",pch=19,cex=1.5,col="green",lwd=2)}
    if (chain[j]=="beta"){  points(x,y,"b",pch=19,cex=1.5,col="blue",lwd=2)}
    
}

imageSave(file=paste0(drive,folder_fig,"individual_dynamics/",ID,"_timecourse.png"))

}
