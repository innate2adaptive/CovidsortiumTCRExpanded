library(ggplot2)
library(reshape)
library(stringr)

estimate_frequencies <- function(df){
  min<-unique(sort(unlist(df),decreasing=FALSE))[2]
  #calculate proportion of each TCR which is zero
  #zeros<-function(x) length(which(data[x,]<min))
  
  pzero<-c()
  j<-1
  for ( j in 1:dim(df)[1]){
    pzero[j]<-length(which(df[j,]<min)) # counts in how many mice a CDR3 is not found
  }
  pzero[pzero == 0]<-0.99 # when present in all mice, to avoid Inf
  i_gt1<-which(pzero<10)
  zeros<-totals[[epitope]]-dim(df)[1]
  # print(zeros)
  #average number of UMI in naive from Michal
  UMI_T<-50000
  mu<--log(pzero[i_gt1]/dim(df)[2])#derive m from proportion of mice which do not have sequence
  f<-10^6*mu/UMI_T
  f_permln<-f/10^6
  names(f_permln)<-rownames(df)
  return(f_permln)
}


gp66eff<-read.csv("data/LCMVdata/E_GP66_freq.csv",header = TRUE,row.names = 1)
gp92eff<-read.csv("data/LCMVdata/E_GP92_freq.csv",header = TRUE,row.names = 1)
np205eff<-read.csv("data/LCMVdata/E_NP205_freq.csv",header = TRUE,row.names = 1)
np396eff<-read.csv("data/LCMVdata/E_NP396_freq.csv",header = TRUE,row.names = 1)

# I want to remove all the ones which have UMI count 0 at day 8 for all mice.

effector_list<-list(GP66=gp66eff, GP92=gp92eff, NP205=np205eff, NP396=np396eff)

## Section B - compare direct and inference methods - This still does not work

input_folder<-"data/LCMVdata/"
dir<-dir("data/LCMVdata/")
dir

#total number of epitope specific TCRs after filtering
totals <-list(GP66=320,GP92=401,NP205=553,NP396=503)

counts_all<-data.frame()
estimates<-data.frame()
effector_timepoints<-data.frame()

filenames<-c("NaiveGP66.csv", "NaiveGP92.csv", "NaiveNP205.csv", "NaiveNP396.csv")

for (f in filenames){
  data1<-read.table(paste0(input_folder,f),sep=",",header = TRUE,row.names = 1)
  
  print(dim(data1))
  print(f)
  
  name<-colnames(data1)[1]
  epitope<-strsplit(name,split="\\.")[[1]][1]
  eff<-effector_list[[epitope]]
  
  ##################################################################
  ## first, day 8 
  ##################################################################
  
  cols<-c(paste0(epitope, ".LCMV8.M1"), paste0(epitope, ".LCMV8.M2"), 
          paste0(epitope, ".LCMV8.M3"), paste0(epitope, ".LCMV8.M4"))
  counts_eff<-rowSums(eff[cols])
  good<-names(counts_eff[counts_eff>0])
  data<-data1[rownames(data1) %in% good,]
  eff_day8<-eff[rownames(eff) %in% good,]
  print("dim day8")
  print(dim(eff_day8))
  colnames(eff_day8)<-lapply(colnames(eff), function(x){
    mylist<-strsplit(x, split="\\.")[[1]][2:3]
    paste(mylist, collapse = ".")
  })
  eff_day8$set<-"day8"
  eff_day8$epitope<-epitope
  not_in_naive<-setdiff(good, rownames(data))

  f_permln<-estimate_frequencies(data)  
  absent<-data.frame(rep(5*10^(-7), length(not_in_naive)))
  rownames(absent)<-not_in_naive
  colnames(absent)<-c("f_permln")
  
  inferred<-rbind(data.frame(f_permln), absent)  

  counts <- data.frame(table(inferred))
  counts$inferred<-as.numeric(as.character(counts$inferred))
  # counts[dim(counts)[1] + 1,]<-c(min(f_permln)/2, zeros)
  counts$log<-log10(counts$inferred)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)

  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/length(good)
  print(counts_agg)  

  counts_agg$epitope<-epitope
  counts_agg$type<-"inferred"
  counts_agg$set<-"day8"
  
  counts_all<-rbind(counts_all, counts_agg)  
  
  # now calculate frequencies directly, rather than from Poisson
  p<-rowSums(data)/10
  names(p)<-rownames(data)
  p1<-data.frame(p)
  colnames(p1)<-c("f_permln")
  direct<-rbind(p1, absent) 
  
  counts <- data.frame(table(direct))
  counts$direct<-as.numeric(as.character(counts$direct))
  
  counts$log<-log10(counts$direct)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/sum(counts_agg$x)
  print(counts_agg)  
  
  counts_agg$epitope<-epitope
  counts_agg$type<-"direct"
  counts_agg$set<-"day8"
  
  counts_all<-rbind(counts_all, counts_agg)
  
  
  colnames(direct)<-c("direct")
  colnames(inferred)<-c("inferred")
  
  both.methods<-merge(direct, inferred, by = 0)
  rownames(both.methods)<-both.methods$Row.names
  both.methods<-both.methods[,2:3]
  both.methods$epitope<-epitope
  both.methods$day<-"day8"
  
  estimates<-rbind(estimates, both.methods)
  
  rm(both.methods)
  rm(p)
  rm(p1)
  rm(data)
  rm(f_permln)
  rm(inferred)
  rm(direct)
  rm(counts)
  rm(counts_agg)
  rm(good)
  rm(not_in_naive)
  rm(absent)
  
  ##################################################################
  ## then day 40
  ##################################################################
  
  eff<-effector_list[[epitope]]
  cols<-c(paste0(epitope, ".LCMV40.M1"), paste0(epitope, ".LCMV40.M2"), 
          paste0(epitope, ".LCMV40.M3"))
  counts_eff<-rowSums(eff[cols])
  good1<-names(counts_eff[counts_eff>0]) # the ones present at day 40
  cols<-c(paste0(epitope, ".LCMV8.M1"), paste0(epitope, ".LCMV8.M2"), 
          paste0(epitope, ".LCMV8.M3"), paste0(epitope, ".LCMV8.M4"))
  counts_eff<-rowSums(eff[cols])
  good2<-names(counts_eff[counts_eff==0]) # the ones not present at day 8
  good<-intersect(good1, good2)
  
  # I want only the ones that come up later
  data<-data1[rownames(data1) %in% good,]
  eff_day40<-eff[rownames(eff) %in% good,]
  print("dim day40")
  print(dim(eff_day40))
  colnames(eff_day40)<-lapply(colnames(eff), function(x){
    mylist<-strsplit(x, split="\\.")[[1]][2:3]
    paste(mylist, collapse = ".")
  })
  eff_day40$set<-"day40"
  eff_day40$epitope<-epitope
  effectors_timepoints<-rbind(effector_timepoints, eff_day8, eff_day40)
  
  not_in_naive<-setdiff(good, rownames(data))
  
  f_permln<-estimate_frequencies(data)  
  absent<-data.frame(rep(5*10^(-7), length(not_in_naive)))
  rownames(absent)<-not_in_naive
  colnames(absent)<-c("f_permln")
  
  inferred<-rbind(data.frame(f_permln), absent)  
  
  counts <- data.frame(table(inferred))
  counts$inferred<-as.numeric(as.character(counts$inferred))

  counts$log<-log10(counts$inferred)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/sum(counts_agg$x)
  print(counts_agg)  
  
  counts_agg$epitope<-epitope
  counts_agg$type<-"inferred"
  counts_agg$set<-"day40"
  
  counts_all<-rbind(counts_all, counts_agg)  
  
  # now calculate frequencies directly, rather than from Poisson
  p<-rowSums(data)/10
  names(p)<-rownames(data)
  p1<-data.frame(p)
  colnames(p1)<-c("f_permln")
  direct<-rbind(p1, absent) 
  
  counts <- data.frame(table(direct))
  counts$direct<-as.numeric(as.character(counts$direct))
  
  counts$log<-log10(counts$direct)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/sum(counts_agg$x)
  print(counts_agg)  
  
  counts_agg$epitope<-epitope
  counts_agg$type<-"direct"
  counts_agg$set<-"day40"
  
  counts_all<-rbind(counts_all, counts_agg)
  
  colnames(direct)<-c("direct")
  colnames(inferred)<-c("inferred")
  
  both.methods<-merge(direct, inferred, by = 0)
  rownames(both.methods)<-both.methods$Row.names
  both.methods<-both.methods[,2:3]
  both.methods$epitope<-epitope
  both.methods$day<-"day40"
  
  estimates<-rbind(estimates, both.methods)
  
  rm(both.methods)
  rm(p)
  rm(p1)
  rm(data)
  rm(f_permln)
  rm(inferred)
  rm(direct)
  rm(counts)
  rm(counts_agg)
  rm(good)
  rm(not_in_naive)
  rm(absent)
  
  rm(data1)
  rm(eff)
}

effectors_timepoints$id<-rownames(effectors_timepoints)

effectors_timepoints_m<-melt(effectors_timepoints)
effectors_timepoints_m$value[effectors_timepoints_m$value == 0] <- 10^(-8)
effectors_timepoints_m$timepoint<-unlist(lapply(as.character(effectors_timepoints_m$variable),
                                         function(x) {
                                           strsplit(x, split="\\.")[[1]][1]
                                         }))
library(dplyr)
effector_timepoints_agg<-effectors_timepoints_m[,c("timepoint", "id", "epitope", "set", "value")] %>% group_by(timepoint, id, epitope, set) %>% summarise(mean_value = mean(value))


effector_timepoints_agg$timepoint<-factor(effector_timepoints_agg$timepoint, 
                                     levels = c("PBS", "LCMV8", "LCMV40"))



p<-ggplot(effector_timepoints_agg) +
  geom_jitter(aes(x=timepoint, y = mean_value, col = set), alpha = 0.5, position = position_dodge(width = 0.1), size = 3) +
  geom_line(aes(x=timepoint, y = mean_value, col = set, group = id), alpha = 0.5) +
  scale_color_manual(values = c(day8 = "darkorange2", day40 = "cyan3")) +
  scale_y_continuous(trans="log10") +
  labs(x = "mice group", y = "mean frequency across mice in group", color = "effector set") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        title=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

svg("output_figures/FigS14A_day8_day40_definition.svg")
print(p)
dev.off()


## then move on to controls

input_c<-"data/LCMVdata/controls/"
ctrl_tot<-1777 # number of sampled sequences

for (i in 1:10){
  filename<-paste0(input_c,"ControlNum", i, ".csv")
  print(filename)
  data<-read.table(filename,sep=",",header = TRUE,row.names = 1)
  print(dim(data))
  epitope<-"control"
  
  f_permln<-estimate_frequencies(data)  
  absent<-data.frame(rep(5*10^(-7), 1777 - dim(data)[1]))
  rownames(absent)<-seq(1:(1777 - dim(data)[1]))
  colnames(absent)<-c("f_permln")
  
  inferred<-rbind(data.frame(f_permln), absent)  
  
  counts <- data.frame(table(inferred))
  counts$inferred<-as.numeric(as.character(counts$inferred))
  counts$log<-log10(counts$inferred)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/sum(counts_agg$x)
  print(counts_agg)  
  
  counts_agg$epitope<-"control"
  counts_agg$type<-"inferred"
  counts_agg$set <-"day8"
  
  counts_all<-rbind(counts_all, counts_agg)
  
  # now calculate frequencies directly, rather than from Poisson
  p<-rowSums(data)/10
  names(p)<-rownames(data)
  p1<-data.frame(p)
  colnames(p1)<-c("f_permln")
  direct<-rbind(p1, absent) 
  
  counts <- data.frame(table(direct))
  counts$direct<-as.numeric(as.character(counts$direct))
  
  counts$log<-log10(counts$direct)
  
  counts$logbin<-unlist(lapply(counts$log, function(x){ifelse(x <= -6, "< 10^-6", 
                                                              ifelse(x <= -5, "10^-6 - 10^-5", 
                                                                     ifelse(x <= -4, "10^-5 - 10^-4", 
                                                                            ifelse(x <= -3, "10^-4 - 10^-3", "> 10^-3"))))}))
  counts_agg<-aggregate(counts$Freq, by=list(logbin = counts$logbin), FUN=sum)
  
  counts_agg$logbin<-factor(counts_agg$logbin, 
                            levels = c("< 10^-6", "10^-6 - 10^-5", "10^-5 - 10^-4", "10^-4 - 10^-3", "> 10^-3"))
  counts_agg$x<-as.numeric(as.character(counts_agg$x))
  counts_agg$x<-counts_agg$x/sum(counts_agg$x)
  print(counts_agg)  
  
  counts_agg$epitope<-epitope
  counts_agg$type<-"direct"
  counts_agg$set<-"day8"
  
  counts_all<-rbind(counts_all, counts_agg)
  
  colnames(direct)<-c("direct")
  colnames(inferred)<-c("inferred")
  
  both.methods<-merge(direct, inferred, by = 0)
  rownames(both.methods)<-both.methods$Row.names
  both.methods<-both.methods[,2:3]
  both.methods$epitope<-"control"
  both.methods$day<-"day8"
  
  estimates<-rbind(estimates, both.methods)
  
  rm(both.methods)
  rm(p)
  rm(p1)
  rm(data)
  rm(f_permln)
  rm(inferred)
  rm(direct)
  rm(counts)
  rm(counts_agg)
  rm(absent)
  
}

counts_all$x<-as.numeric(as.character(counts_all$x))
counts_all_agg<-aggregate(counts_all$x, by=list(epitope=counts_all$epitope, type=counts_all$type, logbin=counts_all$logbin, set=counts_all$set), FUN=mean)


counts_all_agg$type<-factor(counts_all_agg$type, levels = c("inferred", "direct"))
counts_all_agg$set<-factor(counts_all_agg$set, levels = c("all", "day8", "day40"))
counts_all_agg$epitope<-factor(counts_all_agg$epitope, levels = c("GP66", "GP92", "NP205", "NP396", "control"))

p1<-ggplot(counts_all_agg) +
  geom_bar(aes(x = logbin, y = as.numeric(x), fill = set), 
           stat = "identity", alpha = 0.5, position = position_dodge(preserve = "single"), width = 0.7) +
  scale_fill_manual(values = c(day8 = "darkorange2", day40 = "cyan3"))  +
  # scale_color_manual(values = c(inferred = "black", direct = "white")) +
  scale_x_discrete(labels = c(bquote('<'*10^-6),bquote(10^-6*' - '*10^-5), 
                              bquote(10^-5*' - '*10^-4), bquote(10^-4*' - '*10^-3), bquote('>'*10^-3))) +
  labs(x = "TCR frequency", y = "proportion of CDR3s", fill = "") +
  facet_grid(rows = vars(epitope), cols = vars(type), scales = "free") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        title=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)
# ylim(0, 600)


svg("output_figures/Fig5c_sharing_estimate_TCR_LCMV_day8VSday40not8_freq.svg")
print(p1)
dev.off()


p2<-ggplot(estimates, aes(x=direct, y=inferred)) +
  geom_jitter(width = 0.05, height=0.05, alpha=0.2, size=2) +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  labs(x = "UMI-derived CDR3 frequency", y = "Poisson-derived CDR3 frequency") +
  # facet_grid(rows=vars(epitope), cols = vars(day)) +
  geom_abline() +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        title=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=12)) + theme(aspect.ratio = 0.9)

cor.test(estimates$direct, estimates$inferred, method = "pearson")

svg("output_figures/FigS14B_direct_vs_inferred_scatterplot.svg")
print(p2)
dev.off()