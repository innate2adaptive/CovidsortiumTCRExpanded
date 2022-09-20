####Control #####
####libraries####
library(ggplot2)
library(plyr)
library(stringr) 
library(reshape2) # for melt
library(lattice)
library(permute)
library(ggsci)
library(FactoMineR)
library(dplyr)
library(grid)
library(ggrepel)
library(ggpubr)
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))  


#colorblind- palette
cbP_Grey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbP_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbP_Grey_b= c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")






#load data####

#The BPS/LCM 8/LCM 40 and epitope -specific TCR beta data
df.LCMV_sub=readRDS(file = "~/df.LCMV.aa_beta.RDS") #path of the LCMV infected mice
df.LCMV_sub$Class=str_sub(df.LCMV_sub$ClassState,1,3)#add CD4+ CD8+ class 
df.LCMV_sub=subset(df.LCMV_sub,Chain=="B") #only the TCR beta chain 
df.LCMV_sub$Tetramer=str_sub(df.LCMV_sub$Class_Tetramer,5,-5) #Naming Tetramer colomn 
df.LCMV_sub$Tetramer=gsub("GP-6","GP66",df.LCMV_sub$Tetramer)
df.LCMV_sub$classOf_Tet= str_sub(df.LCMV_sub$Class_Tetramer,1,3)
df.LCMV_sub$CDR3_AA_Class= paste(df.LCMV_sub$CDR3_AA,df.LCMV_sub$Class)
  
#Selecting epitope -specific 
tetClass=subset(df.LCMV_sub,classOf_Tet==Class)
length(unique(tetClass$CDR3_AA))#1777

##only naive SP cells from PBS/LCMV8/LCMV40 mice
df.LCMV_subN=subset(df.LCMV_sub,State=="N"&Tissue=="SP")


#Young mice, SP naive_ load data and subset
df.age=readRDS(file = "~/Young_SP_naive.aa_beta.RDS") #path of the SP naive young MICE 
df.age$CDR3_AA_Class=paste(df.age$CDR3_AA,df.age$Class)
length(unique(df.age$CDR3_AA)) #730572


#test how many random control 1777 TCRs are found in the SP naive - PBS/LCMV8/LCMV40 each mouse ##

y2  <- NULL;
for (i in 1:50){
  inds <- sample.int(nrow(df.age),1777, replace=T) #select Randomly 1777 CDR3AA froun young naive TCR
  df.age_s=df.age[inds,]
  dfC= ddply(df.LCMV_subN,.(Mice,Con),function(x){
    #for each mice and condition (PBS, LCMV8,LCMV40) return the overlap CDR3AA with the young mice 
    s1=subset(x,x$CDR3_AA_Class%in%df.age_s$CDR3_AA_Class) 
    CDR3_AA_Class=length(unique(s1$CDR3_AA_Class))
    UMI= sum(s1$sumUMI)
    return(data.frame(CDR3_AA_Class,UMI)) #return the number of CDR3AA and the UMI count in each mice and condition 
  })
  
  y2 <- rbind(y2, dfC)
}

#mean counts over 50 reapets 
df.cont1<- ddply(y2, .(Mice,Con), summarize,
                 UMI_cont = mean(UMI),numCDR3AA_cont=mean(CDR3_AA_Class))

#add the CDR3AA number and UMI counts sum for the epitope -specific ,All and only naive from each mice and condition (PBS,LCMV8,LCMV40)
tetClass_number=ddply(tetClass,.(Mice,Con),function(x)data.frame(numCDR3AA_ES= length(unique(x$CDR3_AA)),UMI_ES=sum(x$sumUMI)))

numall= ddply(df.LCMV_sub,.(Mice,Con),function(x)data.frame(numCDR3AA_ALL= length(unique(x$CDR3_AA)),UMI_ALL=sum(x$sumUMI)))

numallN= ddply(df.LCMV_subN,.(Mice,Con),function(x)data.frame(numCDR3AA_N= length(unique(x$CDR3_AA)),UMI_N=sum(x$sumUMI)))

df2=left_join(tetClass_number,numall,by=c("Mice","Con"))
df3=left_join(df2,df.cont1,by=c("Mice","Con"))
df=left_join(df3,numallN,by=c("Mice","Con"))

#final table in df#
