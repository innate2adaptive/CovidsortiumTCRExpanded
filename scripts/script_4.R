#Annotation of COVID-19 expanded TCRs

library(kernlab)
library(igraph)
library(stringdist)
library(tidyr)
library(dplyr)
library(data.table)
library(pheatmap)

library(circlize)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#on Linux
drive <- "/media/benny/data/"



#on Windows
drive<-"C:/users/Benny/"

#folder for plots
folder_plots<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Figs\\"
#folder for data output
folder_data<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Data_for_paper\\"

#corrected list of epitopes - without TCRs
file1<-"epitope_list.txt"
#matched TCRs
file2<-"exp_AB_merge_2.txt"

epitopes<-read.table(file=paste0(drive, input_data,file1),sep="\t",stringsAsFactors = FALSE,header = TRUE)
matches1<-read.table(file=paste0(drive, input_data,file2),sep="\t",stringsAsFactors = FALSE,header = TRUE)

########################################################################################
#Fig 1c

#early expanded set 
file<-"exp_AB_wide3.RData"
chain<-"alpha"
chain<-"beta"
file_name<-load(paste0(drive,input_data,file))
i_chain<-which(exp_AB_wide3$chain==chain)
i_PCR<-which(exp_AB_wide3$control==FALSE)
length(unique(exp_AB_wide3$junction_aa[i_PCR]))

#remove annotated TCRs which are not in exp_AB_wide3 (may have been removed as MAITs)
both<-match(matches1$junction_aa,unique(exp_AB_wide3$junction_aa[i_PCR]))
matches<-matches1[which(!is.na(both)),]

mat<-matches[,10:40]
i_zero<-which(matches[,10:40]==0,arr.ind = TRUE)
i_one<-which(matches[,10:40]==1,arr.ind = TRUE)

CDR3<-paste0(1:240,rep("_",240),matches$junction_aa)
mat[i_zero]<-"no"
mat[i_one]<-"yes"
rownames(mat)<-CDR3
col_fun1<-c("no" ="gray95","yes"="darkgreen")

#SECTION A - Circular heatmap with annotations
circos.clear()
circos.par(gap.degree = 10)
i<-1
par(mar=c(1,1,1,1))
circos.heatmap(mat[,1:2],cluster=FALSE,col=col_fun1, track.height=0.01, track.margin = c(0.003,0.003),rownames.side = "outside",rownames.cex = 0.2)

for(i in 3:31){
circos.heatmap(mat[,i],cluster=FALSE,col=col_fun1, track.height=0.01, track.margin = c(0.003,0.003))
}

chain<-c("alpha","beta")
col_chain<-c("red", "violet")
names(col_chain)<-chain
circos.heatmap(matches$Chain,col=col_chain,track.height=0.02)

#CD8 CD4
CD<-levels(factor(matches$CD8_CD4))
col_CD<-c("blue","yellow")
names(col_CD)<-CD
circos.heatmap(matches$CD8_CD4,col=col_CD,track.height=0.02)

#HLA
HLA<-levels(factor(matches$MHC))
col_HLA<-rainbow(10)
names(col_HLA)<-HLA
circos.heatmap(matches$MHC,col=col_HLA,track.height=0.02)

#Antigen
Ag<-levels(factor(matches$Epitope.gene))
col_ag<-rainbow(13)
names(col_ag)<-Ag
circos.heatmap(matches$Epitope.gene,col=col_ag,track.height=0.02)

#source of TCRs
source<-levels(factor(matches$Source))
col_source<-rainbow(3)
names(col_source)<-source
circos.heatmap(matches$Source,col=col_source,track.height=0.02)


imageSave(file=paste0(drive,folder_fig,"circus_hm.png"))
circos.clear()

Legd_chain<-Legend(labels = chain ,  legend_gp= gpar(fill=col_chain),title="TCR chain")
Legd_CD<-Legend(labels = CD,  legend_gp= gpar(fill=col_CD),title="CD4 CD8")
Ledg_HLA<-Legend(labels = HLA, legend_gp= gpar(fill= col_HLA), title = "HLA")
Ledg_Ag<-Legend(labels = Ag, legend_gp= gpar(fill= col_ag), title = "Antigen")
Ledg_source<-Legend(labels = source, legend_gp= gpar(fill= col_source), title = "Source")
dev.off()

pd <- packLegend(Legd_chain,Legd_CD, Ledg_HLA, Ledg_Ag, Ledg_source, direction = "horizontal")
draw(pd)

imageSave(file=paste0(drive,folder_fig,"circus_hm_legend.png"))

?circos.heatmap
?colors
?Legend

###################################################################################
#Section B : time course of annotated expanded
match<-as.numeric(sapply(1:dim(exp_AB_wide3)[1],function(x) {exp_AB_wide3$junction_aa[x] %in% unique(matches$junction_aa)}))
i_match<-which(match==1)
exp_AB_ann1<-exp_AB_wide3[i_match,]
PCR_total1<-exp_AB_wide3[i_match,10:18]
#not logged looks nicer !
title<-"time course of annotated expanded TCRs"
plot(-1,-1,xlim=c(0,10),ylim=c(0,(max(PCR_total1,na.rm=TRUE))),main=title,xaxt="none",las=1,xlab="Weeks to PCR+", ylab="TCR abundance (counts/million x 10^3)",yaxt="none",cex.lab=1.5)
#plot(-1,-1,xlim=c(0,10),ylim=c(0,80000),xaxt="none",las=1,xlab="Weeks to PCR+", ylab="Total TCR abundance (counts/million x 10^5)",yaxt="none",cex.lab=1.5)

axis(1,at=c(1:9),label=c(-3:4,15),cex.axis=1.5)
axis(2,at=c(0,1,2,3)*1E3,labels=c(0,1,2,3),cex.axis=1.5,las=1)
i<-1
#c(log2(min/2),8:15)
for ( i in 1 : dim(PCR_total1)[1]){
  i_na<-which(!is.na(PCR_total1[i,]))
  x<-runif(i_na,i_na-(i_na/30),i_na+(i_na/30))
  y<-unlist(PCR_total1[i,i_na])
  points(x,y,pch=19,cex = 1.1,col="blue",type="b")
}

#boxplot((PCR_total1),add=TRUE,col=NA,axes=FALSE, at=1:9,border = "blue",boxlwd=2,pars = list(boxwex = 0.8, staplewex = 1, outwex = 1),notch=F)
imageSave(file=paste0(drive,folder,"timecourse_annotated.png"))
#########
#merge with individual TCR annotation 
exp_AB_merge<-distinct(left_join(exp_AB_ann1[,1:4],matches[,1:8] , by = "junction_aa"))

#Section C : compare V usage of expanded and annotated
#remove those without V annotation
i_V<-which(!is.na(exp_AB_merge$V))
exp_AB_merge_1<-exp_AB_merge[i_V,]
i<-1
V_match<-sapply(1:dim(exp_AB_merge_1)[1], function(i) {grepl(exp_AB_merge_1$v_call[i],exp_AB_merge_1$V[i])},simplify=TRUE)
V_match_l<-sum(V_match)



#as a control do random permutation
true<-c()
for (k in 1 :10){
  i_s<-sample(1:dim(exp_AB_merge)[1],replace=TRUE)
  exp_AB_merge_s<-exp_AB_merge_1[i_s,]
  V_match<-sapply(1:dim(exp_AB_merge_1)[1], function(i) {grepl(exp_AB_merge_s$v_call[i],exp_AB_merge_1$V[i])},simplify=TRUE)
  true[k]<-sum(V_match)
  
}
V_match_ctrl<-mean(true)

#Fisher
a<-c(V_match_l,length(V_match)-V_match_l)
b<-c(V_match_ctrl,length(V_match)-V_match_ctrl)
m<-rbind(a,b)
fisher.test(m)

par(lwd=2,mar=c(5,5,2,5))
barplot(c(V_match_l,V_match_ctrl), width = c(0.5,0.5), beside = TRUE, main = paste(" V gene"),density=c(20,0),names.arg = c("Expanded", "Control"),  col=c("blue","black"),lwd=2, ylab="Matched V genes",cex.lab=1.5, cex.names=1.5,xlim = c(0,2))
imageSave(file=paste0(drive,folder,"V_usage.png"))


#CI<-signif(fisher.test(m)[[2]],2)
#ODD<-signif(fisher.test(m)[[3]],2)
#p<-signif(fisher.test(m)[[1]],3)
#par(lwd=2)
#barplot(t(m), width = c(0.5,0.5), beside = FALSE, main = paste(ODD,"\n",CI[1],CI[2],"\n",p),density=c(20,0),names.arg = c("Observed","Control"),  col=c("blue","black"),lwd=2, yaxt="n" ,ylab="Matched V gene usage",cex.lab=1.5, legend.text = c("Match", "No match"),args.legend = list(x = 0.3, y = -10),cex.names=1.5,xlim = c(0,2))
#axis(2, at = c(0,50,100,150,200,250), labels=c(0,50,100,150,200,250),cex.axis=1.5,las=1)
#imageSave(file=paste0(drive,folder,"V_usage.png"))

#Section D : compare HLA of annotated and expanded individuals
HLA<-read.table(file=paste0(drive,input_data,"summary_with_HLA.txt"),header = TRUE, sep="\t")

#convert thw annoated into long form
matches_long_1<-pivot_longer(matches,cols=c(10:40),names_to="ID",values_to = "found",names_prefix = "X")
i_found<-which(matches_long_1$found==1)
matches_long_2<-distinct(matches_long_1[i_found,])
matches_long_2$ID<-as.numeric(matches_long_2$ID)
exp_AB_HLA<-distinct(left_join(matches_long_2,HLA,by = "ID", keep = TRUE))

#remove rows without a known MHC for annotated TCR
i_HLA<-which(!is.na(exp_AB_HLA$MHC))
exp_AB_HLA1<-exp_AB_HLA[i_HLA,]
match<-c()
i<-1
for (i in 1:dim(exp_AB_HLA1)[1]){
MHC<-exp_AB_HLA1$MHC[i]
MHC_1<- gsub ("\\*","\\\\*",MHC)
match[i]<-length(grep(MHC_1,exp_AB_HLA1[i,23:34]))
}

m<-length(which(match>0))       

nm<-length(which(match==0))     
m+nm
#controls 

m_c<-c()
for (k in 1:10){
  
i_s<-sample(1:dim(exp_AB_HLA1)[1],replace=TRUE)
exp_AB_HLA2<-exp_AB_HLA1[i_s,]
i<-1
for (i in 1:dim(exp_AB_HLA1)[1]){
MHC<-exp_AB_HLA2$MHC[i]
MHC_1<- gsub ("\\*","\\\\*",MHC)
match[i]<-length(grep(MHC_1,exp_AB_HLA1[i,23:34]))
}
m_c[k]<-length(which(match>0))       
}
m_c<-mean(m_c)
nm_c<-length(i_HLA)-mean(m_c)
a<-c(m,nm)
b<-c(m_c,nm_c)
c<-rbind(a,b)
fisher.test(c)

par(lwd=2,mar=c(5,5,2,5))
barplot(c(m,m_c), width = c(0.5,0.5), beside = TRUE, main = paste(" HLA restriction"),density=c(20,0),names.arg = c("Expanded", "Control"),  col=c("blue","black"),lwd=2, ylab="Matched HLA",cex.lab=1.5, cex.names=1.5,xlim = c(0,2))
imageSave(file=paste0(drive,folder,"HLA_usage.png"))

CI<-signif(fisher.test(c)[[2]],2)
ODD<-signif(fisher.test(c)[[3]],2)
p<-signif(fisher.test(c)[[1]],3)
par(lwd=2)
barplot(t(c), width = c(0.5,0.5), beside = FALSE, main = paste(ODD,"\n",CI[1],CI[2],"\n",p),density=c(20,0),names.arg = c("Observed","Control"),  col=c("blue","black"),lwd=2, yaxt="n" ,ylab="Matched HLA gene usage",cex.lab=1.5, legend.text = c("Match", "No match"),args.legend = list(x = 0.3, y = -10),cex.names=1.5,xlim = c(0,2))
axis(2, at = c(0,50,100,150,200,250), labels=c(0,50,100,150,200,250),cex.axis=1.5,las=1)
imageSave(file=paste0(drive,folder,"HLA_match.png"))
