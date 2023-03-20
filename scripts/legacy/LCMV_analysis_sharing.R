#
library(mgcv)
library(tidyr)
library(dplyr)
library(zoo)
library(kernlab)
library(igraph)
library(stringdist)
library(MASS)
library(viridis)
library(KernSmooth)

input_folder<-"data/LCMVdata/"
dir<-dir("data/LCMVdata/")
dir
#total number of epitope specific TCRs afeter filtering
totals <-c(738,782,600,610)

i<-1
k<-1
for (i in c(1:4)){
print(dir[i+4])
data<-read.table(paste0(input_folder,dir[i+4]),sep=",",header = TRUE,row.names = 1)
data_E<-read.table(paste0(input_folder,dir[i]),sep=",",header = TRUE,row.names = 1)

name<-colnames(data)[1]
epitope<-strsplit(name,split="\\.")[[1]][1]
 #add columns and look for differences between columns
colsums<-colSums(data)
colsums_E<-colSums((data_E))
colsums_E
colsums

x<-c(3,3,3,2,2,2,2,1,1,1)
anova<-anova(lm(colsums~x))
p<-anova[[5]][1]
anova

lm(colsums_E~x)
anova<-anova(lm(colsums_E~x))
p<-anova[[5]][1]
anova
print(wilcox.test(colsums_E[4:7],colsums))


x_jit<-rnorm(10,x,0.1)

par(mar=c(5,8,2,5))


plot(x_jit+0.5,colsums_E,pch=19,xlim=c(0,4),main=epitope,ylab="TCR frequency",xlab="",xaxt="n",yaxt="n")
points(x_jit,colsums,pch=19,col="blue")
axis(1,at=c(1,2,3),labels=c("PBS","Day 8","Day 40"),cex.axis=1.2)
axis(2,las=2,cex.axis=1)
#?text
#text(0.2,0.009,paste0("P= ",round(p,2)))
}
# folder

min<-unique(sort(unlist(data),decreasing=FALSE))[2]
#calculate proportion of each TCR which is zero
#zeros<-function(x) length(which(data[x,]<min))

pzero<-c()
j<-1
for ( j in 1:dim(data)[1]){
  pzero[j]<-length(which(data[j,]<min))
}
i_gt1<-which(pzero<10)
zeros<-totals[k]-dim(data)[1]
#average number of UMI in naive from Michal
UMI_T<-50000
mu<--log(pzero[i_gt1]/dim(data)[2])

p_poiss<-mu*1E6/UMI_T

#plot freqeuncy diagram
breaks<-c(1,10,100,1000)
density<-hist(p_poiss,breaks,freq=FALSE)$counts
density

par(mar=c(3,15,5,2))
mp<-barplot(c(zeros,density), col = "blue", las=2,cex.axis=2.5,main=epitope)

imageSave(file=paste0(drive,folder,epitope,"_TCR_poiss.png"))


#######################################################################
#######################################################################
#Calcualte freqeuncies directly, rather from Poissson.

p<-rowSums(data)*1E6/10

#mu<--log(pzero[i_gt1]/dim(data)[2])*1E6/UMI_T
#mu
#plot freqeuncy diagram
breaks<-c( 0,1,10,100,1000)
density<-hist(p,breaks,freq=FALSE)$counts
density
par(mar=c(3,15,5,2))
mp<-barplot(c(zeros+density[1],density[2:4]), col = "blue", las=2,cex.axis=2.5,main=epitope)
#############################################################
#plot controls
mp<-barplot(c(1777,12,3,0), col = "blue", las=2,cex.axis=2.5,main="control")
imageSave(file=paste0(drive,folder,"control_TCR_freq.png"))
##############################################################
#axis(1,at=mp,labels=c(expression(10^-6~),2,3))
#?expression
imageSave(file=paste0(drive,folder,epitope,"_TCR_freq.png"))
k<-k+1
cat(k,"\n")
}
?barplot
