library(ggplot2)
library(reshape)
library(stringr)

gp66eff<-read.csv("data/LCMVdata/E_GP66_UMI.csv",header = TRUE,row.names = 1)
gp92eff<-read.csv("data/LCMVdata/E_GP92_UMI.csv",header = TRUE,row.names = 1)
np205eff<-read.csv("data/LCMVdata/E_NP205_UMI.csv",header = TRUE,row.names = 1)
np396eff<-read.csv("data/LCMVdata/E_NP396_UMI.csv",header = TRUE,row.names = 1)

# I want to remove all the ones which have UMI count 0 at day 8 for all mice.

effector_list<-list(GP66=gp66eff, GP92=gp92eff, NP205=np205eff, NP396=np396eff)

## Section B - compare direct and inference methods - This still does not work

input_folder<-"data/LCMVdata/"
dir<-dir("data/LCMVdata/")
dir
#total number of epitope specific TCRs after filtering
totals <-c(738,782,600,610)
counts_all<-data.frame()

filenames<-c("NaiveGP66.csv", "NaiveGP92.csv", "NaiveNP205.csv", "NaiveNP396.csv")

i<-1
k<-1
for (f in filenames){
  data1<-read.table(paste0(input_folder,f),sep=",",header = TRUE,row.names = 1)
  # data_E<-read.table(paste0(input_folder,dir[i]),sep=",",header = TRUE,row.names = 1)
  
  print(dim(data1))
  
  name<-colnames(data1)[1]
  epitope<-strsplit(name,split="\\.")[[1]][1]
  
  eff<-effector_list[[epitope]]
  cols<-c(paste0(epitope, ".LCMV8.M1"), paste0(epitope, ".LCMV8.M2"), 
          paste0(epitope, ".LCMV8.M3"), paste0(epitope, ".LCMV8.M4"))
  counts_eff<-rowSums(eff[cols])
  good<-names(counts_eff[counts_eff>0])
  data<-data1[rownames(data1) %in% good,]
  
  print(epitope)
  print("proportions of effectors kept")
  print(length(good)/dim(eff)[1])
  print("proportions of naive kept")
  print(dim(data)[1]/dim(data1)[1])
  
  min<-unique(sort(unlist(data),decreasing=FALSE))[2]
  #calculate proportion of each TCR which is zero
  #zeros<-function(x) length(which(data[x,]<min))
  
  pzero<-c()
  j<-1
  for ( j in 1:dim(data)[1]){
    pzero[j]<-length(which(data[j,]<min)) # counts in how many mice a CDR3 is not found
  }
  i_gt1<-which(pzero<10)
  zeros<-totals[k]-dim(data)[1]
  # print(zeros)
  #average number of UMI in naive from Michal
  UMI_T<-50000
  mu<--log(pzero[i_gt1]/dim(data)[2])#derive m from proportion of mice which do not have sequence
  f<-10^6*mu/UMI_T
  f_permln<-f/10^6
  
  counts <- data.frame(table(f_permln))
  counts$f_permln<-as.numeric(as.character(counts$f_permln))
  # counts[dim(counts)[1] + 1,]<-c(min(f_permln)/2, zeros)
  counts$log<-log10(counts$f_permln)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  if ("< 10^-6" %in% counts_agg$logbin){
    all_small<-as.numeric(counts_agg[counts_agg$logbin == "< 10^-6",]$x) + zeros
  }else{
    all_small<-zeros
  }
  
  print(all_small)
  counts_agg<-counts_agg[counts_agg$logbin != "< 10^-6",]
  counts_agg<-rbind(counts_agg, c("< 10^-6", all_small))
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/totals[k]
  print(counts_agg)
  
  counts_agg$epitope<-epitope
  counts_agg$type<-"inferred"
  
  counts_all<-rbind(counts_all, counts_agg)
  rm(counts_agg)
  rm(mu)
  rm(f)
  rm(f_permln)
  
  #######################################################################
  #######################################################################
  #Calculate frequencies directly, rather from Poisson.
  
  p<-rowSums(data)/10
  
  counts <- data.frame(table(p))
  counts$p<-as.numeric(as.character(counts$p))
  # counts[dim(counts)[1] + 1,]<-c(min(p)/10, zeros)
  counts$log<-log10(counts$p)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4",
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  counts_agg[counts_agg$logbin == "< 10^-6",]$x <- counts_agg[counts_agg$logbin == "< 10^-6",]$x + zeros
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-counts_agg$x/totals[k]
  counts_agg$epitope<-epitope
  counts_agg$type<-"direct"
  
  counts_all<-rbind(counts_all, counts_agg)
  
  # rm(counts_agg)
  
  k<-k+1
  cat(k,"\n")
  rm(data)
  rm(p)
}

# for controls we take average of 10 controls

input_c<-"data/LCMVdata/controls/"
ctrl_tot<-1777

for (i in 1:10){
  filename<-paste0(input_c,"ControlNum", i, ".csv")
  print(filename)
  data<-read.table(filename,sep=",",header = TRUE,row.names = 1)
  print(dim(data))
  epitope<-"control"
  
  min<-unique(sort(unlist(data),decreasing=FALSE))[2]
  #calculate proportion of each TCR which is zero
  #zeros<-function(x) length(which(data[x,]<min))
  
  pzero<-c()
  j<-1
  for ( j in 1:dim(data)[1]){
    pzero[j]<-length(which(data[j,]<min)) # counts in how many mice a CDR3 is not found
  }
  i_gt1<-which(pzero<10)
  zeros<-ctrl_tot-dim(data)[1]
  print(zeros)
  # print(zeros)
  #average number of UMI in naive from Michal
  UMI_T<-50000
  mu<--log(pzero[i_gt1]/dim(data)[2])#derive m from proportion of mice which do not have sequence
  f<-10^6*mu/UMI_T
  f_permln<-f/10^6
  
  counts <- data.frame(table(f_permln))
  counts$f_permln<-as.numeric(as.character(counts$f_permln))
  # counts[dim(counts)[1] + 1,]<-c(min(f_permln)/2, zeros)
  counts$log<-log10(counts$f_permln)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  if ("< 10^-6" %in% counts_agg$logbin){
    all_small<-as.numeric(counts_agg[counts_agg$logbin == "< 10^-6",]$x) + zeros
  }else{
    all_small<-zeros
  }
  
  print(all_small)
  counts_agg<-counts_agg[counts_agg$logbin != "< 10^-6",]
  counts_agg<-rbind(counts_agg, c("< 10^-6", all_small))
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/ctrl_tot
  print(counts_agg)
  
  counts_agg$epitope<-epitope
  counts_agg$type<-"inferred"
  
  counts_all<-rbind(counts_all, counts_agg)
  rm(counts_agg)
  rm(mu)
  rm(f)
  rm(f_permln)
  
  #######################################################################
  #######################################################################
  #Calculate frequencies directly, rather from Poisson.
  
  p<-rowSums(data)/10
  
  counts <- data.frame(table(p))
  counts$p<-as.numeric(as.character(counts$p))
  # counts[dim(counts)[1] + 1,]<-c(min(p)/10, zeros)
  counts$log<-log10(counts$p)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4",
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  counts_agg[counts_agg$logbin == "< 10^-6",]$x <- counts_agg[counts_agg$logbin == "< 10^-6",]$x + zeros
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-counts_agg$x/ctrl_tot
  counts_agg$epitope<-epitope
  counts_agg$type<-"direct"
  
  counts_all<-rbind(counts_all, counts_agg)
  
  # rm(counts_agg)
  
  k<-k+1
  cat(k,"\n")
  rm(data)
  rm(p)
  
}

counts_all$x<-as.numeric(as.character(counts_all$x))
counts_all_agg<-aggregate(counts_all$x, by=list(epitope=counts_all$epitope, type=counts_all$type, logbin=counts_all$logbin), FUN=mean)


counts_all_agg$type<-factor(counts_all_agg$type, levels = c("inferred", "direct"))
counts_all_agg$epitope<-factor(counts_all_agg$epitope, levels = c("GP66", "GP92", "NP205", "NP396", "control"))

p1<-ggplot(counts_all_agg) +
  geom_bar(aes(x = logbin, y = as.numeric(x), fill = type), stat = "identity", alpha = 0.7, color = "black") +
  scale_fill_manual(values = c(inferred = "burlywood", direct = "azure2")) +
  # scale_color_manual(values = c(inferred = "black", direct = "white")) +
  scale_x_discrete(labels = c(bquote('<'*10^-6),bquote(10^-6*' - '*10^-5), 
                              bquote(10^-5*' - '*10^-4), bquote(10^-4*' - '*10^-3), bquote('>'*10^-3))) +
  labs(x = "TCR frequency", y = "number of CDR3s", fill = "") +
  facet_grid(rows = vars(epitope), cols = vars(type)) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 1.3)
# ylim(0, 600)

svg("output_figures/Fig5c_sharing_estimate_TCR_LCMV_day8.svg")
print(p1)
dev.off()
