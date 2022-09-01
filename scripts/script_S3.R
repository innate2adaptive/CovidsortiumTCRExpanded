#Annotation of COVID-19 expanded TCRs
library(kernlab)
library(igraph)
library(stringdist)
library(tidyr)
library(dplyr)
library(data.table)
library(pheatmap)

options(timeout = max(1000, getOption("timeout"))) # increasing timeout might be necessary to load the files

myURL<-"https://www.dropbox.com/s/a7ymcecnpomge2e/all_A_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)
#plot number of reads per sample
abundance<-all_A_long$total*all_A_long$proportion/1E6
sample_ID<-paste(all_A_long$week_PCR,all_A_long$ID,sep="_")
sample_total<-aggregate(abundance,by=list(sample_ID),sum )

rm(all_A_long)
#betas
myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)
#plot number of reads per sample
abundance_B<-all_B_long$total*all_B_long$proportion/1E6
sample_ID_B<-paste(all_B_long$week_PCR,all_B_long$ID,sep="_")
sample_total_B<-aggregate(abundance_B,by=list(sample_ID_B),sum )
rm(all_B_long)

x<-c(rnorm(dim(sample_total)[1],1,0.1),rnorm(dim(sample_total_B)[1],2,0.1))
y<-c(log10(sample_total$x),log10(sample_total_B$x))

par(yaxt="s")
plot(x,y,xlim=c(0,3),ylim=c(3,6),xaxt="n", pch=19,yaxt="n",xlab="",ylab="",cex=0.3)
boxplot(list(log10(sample_total$x),log10(sample_total_B$x)),add=TRUE,col=NA,cex.names = 1.5, names=NA ,yaxt="n")
axis(2,at = c(3,4,5),labels=c("1000","10000","100000"),cex.axis=1.5,las=2)
axis(1, at=c(1,2), labels= c("alpha","beta"),cex.axis=1.5)
imageSave(file=paste0(drive,folder_fig,"alpha&beta_total.png"))

rownames(sample_total)<-sample_total$Group.1
rownames(sample_total_B)<-sample_total_B$Group.1

x_1<-log10(sample_total[intersect(sample_total$Group.1,sample_total_B$Group.1),2])
y_1<-log10(sample_total_B[intersect(sample_total$Group.1,sample_total_B$Group.1),2])

10^median(x_1)
10^median(y_1)

par(mar=c(5,4,4,2))
par(mar=c(5,8,4,2))
plot(x_1,y_1,pch=18,xlab = "TCR alpha",ylab = "TCR beta",ylim=c(4,5.4),xlim=c(4,5.4),xaxt="n",yaxt="n",cex.lab=2,main="0.94")
axis(1,at=c(4,5),c("10,000","100,000"),cex.axis=1.5)
axis(2,at=c(4,5),c("10,000","100,000"),cex.axis=1.5,las=2)
imageSave(file=paste0(drive,folder_fig,"alpha_v_beta.png"))
cor.test(x_1,y_1, method = "spearman")
