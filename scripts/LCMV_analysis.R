library(ggplot2)
library(reshape)
library(stringr)

gp66naive<-read.csv("data/LCMVdata/NaiveGP66_UMI.csv")
gp66eff<-read.csv("data/LCMVdata/E_GP66_UMI.csv")

gp92naive<-read.csv("data/LCMVdata/NaiveGP92_UMI.csv")
gp92eff<-read.csv("data/LCMVdata/E_GP92_UMI.csv")

np205naive<-read.csv("data/LCMVdata/NaiveNP205_UMI.csv")
np205eff<-read.csv("data/LCMVdata/E_NP205_UMI.csv")

np396naive<-read.csv("data/LCMVdata/NaiveNP396_UMI.csv")
np396eff<-read.csv("data/LCMVdata/E_NP396_UMI.csv")

gp66naive_c<-colSums(gp66naive[,2:11])
gp66eff_c<-colSums(gp66eff[,2:11])

gp92naive_c<-colSums(gp92naive[,2:11])
gp92eff_c<-colSums(gp92eff[,2:11])

np205naive_c<-colSums(np205naive[,2:11])
np205eff_c<-colSums(np205eff[,2:11])

np396naive_c<-colSums(np396naive[,2:11])
np396eff_c<-colSums(np396eff[,2:11])

## Section A - total counts (Fig 5B)

counts_gp66<-data.frame(t(rbind(gp66naive_c,gp66eff_c)))
colnames(counts_gp66)<-c("naive", "effector")
counts_gp66$id<-rownames(counts_gp66)
counts_gp66[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_gp66$id, "[.]", 3)

g<-counts_gp66[(counts_gp66$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_gp66$naive))

counts_gp66_m<-melt(counts_gp66)

counts_gp92<-data.frame(t(rbind(gp92naive_c,gp92eff_c)))
colnames(counts_gp92)<-c("naive", "effector")
counts_gp92$id<-rownames(counts_gp92)
counts_gp92[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_gp92$id, "[.]", 3)

g<-counts_gp92[(counts_gp92$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_gp92$naive))

counts_gp92_m<-melt(counts_gp92)

counts_np205<-data.frame(t(rbind(np205naive_c,np205eff_c)))
colnames(counts_np205)<-c("naive", "effector")
counts_np205$id<-rownames(counts_np205)
counts_np205[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_np205$id, "[.]", 3)

g<-counts_np205[(counts_np205$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_np205$naive))

counts_np205_m<-melt(counts_np205)

counts_np396<-data.frame(t(rbind(np396naive_c,np396eff_c)))
colnames(counts_np396)<-c("naive", "effector")
counts_np396$id<-rownames(counts_np396)
counts_np396[c("epitope","LCMV", "mouse")]<-str_split_fixed(counts_np396$id, "[.]", 3)

g<-counts_np396[(counts_np396$LCMV == "LCMV8"),]
print(wilcox.test(g$effector,counts_np396$naive))

counts_np396_m<-melt(counts_np396)

counts_m<-rbind(counts_gp66_m, counts_gp92_m, counts_np205_m, counts_np396_m)


counts_m$LCMV<-factor(counts_m$LCMV, levels = c("PBS", "LCMV8", "LCMV40"))

p<-ggplot(counts_m) + 
  geom_point(aes(x = LCMV, y = value, color = variable), 
              position = position_dodge(width = 0.5), size = 5, alpha = 0.4) +
  labs(y = "TCR abundance", x = "", col = "") +
  scale_x_discrete(labels = c("PBS", "day 8", "day 40")) +
  scale_colour_manual(values = c(effector = "darkseagreen", naive = "darkmagenta")) +
  facet_wrap(vars(epitope), ncol=2) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=14),
        strip.text = element_text(size=16))

svg("output_figures/Fig5b_abundance_TCR_LCMV.svg")
print(p)
dev.off()

## Section B - compare direct and inference methods - This still does not work

input_folder<-"data/LCMVdata/"
dir<-dir("data/LCMVdata/")
dir
#total number of epitope specific TCRs afeter filtering
totals <-c(738,782,600,610)
counts_all<-data.frame()

i<-1
k<-1
for (i in c(1,3,5,7)){
  data<-read.table(paste0(input_folder,dir[i+8]),sep=",",header = TRUE,row.names = 1)
  data_E<-read.table(paste0(input_folder,dir[i]),sep=",",header = TRUE,row.names = 1)
    
  print(dim(data_E))
  
  name<-colnames(data)[1]
  epitope<-strsplit(name,split="\\.")[[1]][1]

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
  mu<--log(pzero[i_gt1]/dim(data)[2])#derive m from proportion of mice which do not have sequence
  f<-10^6*mu/UMI_T
  f_permln<-f/10^6
  
  counts <- data.frame(table(f_permln))
  counts$f_permln<-as.numeric(as.character(counts$f_permln))
  counts[dim(counts)[1] + 1,]<-c(min(f_permln)/10, zeros)
  counts$log<-log10(counts$f_permln)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4",
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts$logbin<-factor(counts$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  
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
  counts[dim(counts)[1] + 1,]<-c(min(p)/10, zeros)
  counts$log<-log10(counts$p)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4",
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts$logbin<-factor(counts$logbin, 
                        levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  
  counts_agg$epitope<-epitope
  counts_agg$type<-"direct"
  
  counts_all<-rbind(counts_all, counts_agg)
  
  rm(counts_agg)
  
  #############################################################
  # #plot controls
  # mp<-barplot(c(1777,12,3,0), col = "blue", las=2,cex.axis=2.5,main="control")
  ##############################################################

  k<-k+1
  cat(k,"\n")
  rm(data)
  rm(p)
}

counts_all$type<-factor(counts_all$type, levels = c("inferred", "direct"))

p1<-ggplot(counts_all) +
  geom_bar(aes(x = logbin, y = x, fill = type), stat = "identity", alpha = 0.7, color = "black") +
  scale_fill_manual(values = c(inferred = "burlywood", direct = "azure2")) +
  # scale_color_manual(values = c(inferred = "black", direct = "white")) +
  scale_x_discrete(labels = c(bquote('<'*10^-6),bquote(10^-6*' - '*10^-5), 
                              bquote(10^-5*' - '*10^-4), bquote(10^-4*' - '*10^-3), bquote('>'*10^-3))) +
  labs(x = "", y = "number of CDR3s", fill = "") +
  facet_grid(rows = vars(epitope), cols = vars(type)) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) 

print(p1)
