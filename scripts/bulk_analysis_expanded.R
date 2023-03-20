# library(tabula)

load("data/output_data/exp_AB_wide3.RData")

idcols <-c("decombinator_id", "v_call", "junction_aa", "j_call", "chain",
           "ID", "control")
variablecol1<-c("total_-3", "total_-2", "total_-1", "total_0", "total_1", "total_2", "total_3",
               "total_4", "total_14")
variablecol2<-c("proportion_-3", "proportion_-2", "proportion_-1", "proportion_0", "proportion_1", 
                "proportion_2",  "proportion_3", "proportion_4", "proportion_14")

exp_AB_long1<-reshape(exp_AB_wide3[c(idcols, variablecol1)], direction = "long", 
        varying = variablecol1, timevar = "week_PCR", 
        times = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "14"),
        v.names = "total", sep = "_")
exp_AB_long2<-reshape(exp_AB_wide3[c(idcols, variablecol2)], direction = "long", 
                      varying = variablecol2, timevar = "week_PCR", 
                      times = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "14"),
                      v.names = "proportion", sep = "_")

exp_AB_long<-merge(exp_AB_long1, exp_AB_long2, by=c(idcols, "week_PCR"))
exp_AB_long[is.na(exp_AB_long)]<-0
exp_AB_long0<-exp_AB_long[exp_AB_long$total > 0,]

exp_AB_long0$week_PCR<-lapply(exp_AB_long0$week_PCR, function(x) {gsub("\\.", "-", x)})
exp_AB_long0$mycounts<-as.numeric(exp_AB_long0$total) * (as.numeric(exp_AB_long0$proportion)/10^6)
exp_AB_long0$counts<-as.integer(round(exp_AB_long0$mycounts, 0))

## note that "total" is not total number of expanded, but rather total number in sample
library(tidyr)
library(dplyr)
library(testit)
library(vegan)
library(DescTools)

exp_AB_long0$total_exp<-0

for (id in unique(exp_AB_long0$ID)){
  for (chain in c("alpha", "beta")){
    ss<-exp_AB_long0[(exp_AB_long0$ID == id) & (exp_AB_long0$chain == chain),]
    for (week in unique(ss$week_PCR)){
      total<-sum(exp_AB_long0[(exp_AB_long0$ID == id) & 
                                  (exp_AB_long0$chain == chain) & 
                                  (exp_AB_long0$week_PCR == week), "counts"])
      exp_AB_long0[(exp_AB_long0$ID == id) & 
                     (exp_AB_long0$chain == chain) & 
                     (exp_AB_long0$week_PCR == week),]$total_exp<-total
    }
  }
}

exp_AB_long01<-exp_AB_long0[(exp_AB_long0$total_exp > 100) & (exp_AB_long0$counts > 0) ,]

exp_A_long<-exp_AB_long01[exp_AB_long01$chain == "alpha",]
exp_B_long<-exp_AB_long01[exp_AB_long01$chain == "beta",]

#min number to subsample to
min_depth<-min(exp_A_long$total_exp) # total contains the total number of sequences from each week
subsampleSize_a<-round(min_depth/10, 0)*10

extract_diversity_metrics<-function(x, expanded_x){
  total<-sum(x$counts)
  richness_cdr3<-length(unique(x$junction_aa))/total
  richness_dcb<-length(unique(x$decombinator_id))/total
  # y<-x[!duplicated(x[,c("junction_aa", "v_call", "j_call")]),]
  x$tcr<-paste(x$v_call, x$junction_aa, x$j_call, sep="-")
  expanded_x$tcr<-paste(expanded_x$v_call, expanded_x$junction_aa, expanded_x$j_call, sep="-")
  richness_tcr<-length(unique(x$tcr))/total
  
  shannon_dcb<-diversity(x$counts, index = "shannon", MARGIN = 1, base = exp(1))
  # to calculate shannon on actual aas, I first need to collapse on actual cdr3 aa counts
  cdr3_unique<-table(expanded_x$junction_aa)/total
  tcr_unique<-table(expanded_x$tcr)/total
  shannon_cdr3<-diversity(cdr3_unique, index = "shannon", MARGIN = 1, base = exp(1))
  shannon_tcr<-diversity(tcr_unique, index = "shannon", MARGIN = 1, base = exp(1))
  gini_cdr3<-Gini(cdr3_unique)
  gini_tcr<-Gini(tcr_unique)
  gini_dcb<-Gini(as.vector(x$counts/total))
  renyi_orders<-c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf)
  renyi_cdr3<-renyi(cdr3_unique, scales = renyi_orders)
  names(renyi_cdr3)<-c("renyi_0_cdr3", "renyi_0.25_cdr3", "renyi_0.5_cdr3", "renyi_1_cdr3", "renyi_2_cdr3",
                       "renyi_4_cdr3", "renyi_8_cdr3", "renyi_16_cdr3", "renyi_32_cdr3", "renyi_64_cdr3", "renyi_Inf_cdr3")
  renyi_dcb<-renyi(as.vector(x$counts)/total, scales = renyi_orders)
  names(renyi_dcb)<-c("renyi_0_dcb", "renyi_0.25_dcb", "renyi_0.5_dcb", "renyi_1_dcb", "renyi_2_dcb",
                       "renyi_4_dcb", "renyi_8_dcb", "renyi_16_dcb", "renyi_32_dcb", "renyi_64_dcb", "renyi_Inf_dcb")
  renyi_tcr<-renyi(tcr_unique, scales = renyi_orders)
  names(renyi_tcr)<-c("renyi_0_tcr", "renyi_0.25_tcr", "renyi_0.5_tcr", "renyi_1_tcr", "renyi_2_tcr",
                       "renyi_4_tcr", "renyi_8_tcr", "renyi_16_tcr", "renyi_32_tcr", "renyi_64_tcr", "renyi_Inf_tcr")
  
  # I should also do this at the tcr level, so I can encompass v and j
  
  my_results<-c("richness_dcb"=richness_dcb, "richness_cdr3"=richness_cdr3, "richness_tcr"=richness_tcr, 
                "shannon_dcb"=shannon_dcb, "shannon_cdr3"=shannon_cdr3, "shannon_tcr"=shannon_tcr,
                "gini_dcb"=gini_dcb, "gini_cdr3"=gini_cdr3, "gini_tcr"=gini_tcr, renyi_cdr3, 
                renyi_tcr, renyi_dcb)
  
  return(my_results)
  
}

results_A<-data.frame()
results_ss_summary_A<-data.frame()

for (id in unique(exp_A_long$ID)){
  print(paste("ID:", id))
  thisID<-exp_A_long[exp_A_long$ID == id,]
  sum(thisID$counts)
  for (week in unique(thisID$week_PCR)){
    print(paste("week:", week))
    thisWeek<-thisID[(thisID$week_PCR == week) & (thisID$counts > 0),]
    print(dim(thisWeek))
    control<-unique(thisWeek$control)
    total<-unique(thisWeek$total_exp)
 
    expanded<-thisWeek %>% uncount(counts)
    assert("Assert expected number of rows in table", total == dim(expanded)[1])
    
    myresults<-extract_diversity_metrics(thisWeek, expanded)
    
    week_results<-c("id"=id, "control" = control, "week"=week, "total"=total, 
                    myresults)
    results_A<-rbind(results_A, t(week_results))
    print("Total results computed")
    
    ss_results_all<-data.frame()
    print("Iterating subsamples")
    t0<-Sys.time()
    
    for (i in 1:100){
      # t<-Sys.time()
      subsampled_expanded<-data.frame(expanded[sample(nrow(expanded), subsampleSize_a), ])
      total_ss<-subsampleSize_a
      
      subsampled<- subsampled_expanded %>% group_by(across(everything())) %>% count()
      subsampled$counts<-subsampled$n
      assert("Collapsing table correctly", dim(subsampled)[1] == dim(subsampled_expanded[!duplicated(subsampled_expanded),])[1])
      
      myresults_ss<-extract_diversity_metrics(subsampled, subsampled_expanded)
      ss_results_all<-rbind(ss_results_all, t(myresults_ss))
      
      # print(Sys.time()-t)
    }
    
    print("Subsampling done")
    print(Sys.time()-t0)
    
    ss_results_summ<-colMeans(ss_results_all)
    
    week_results_ss<-c("id"=id, "control" = control, "week"=week, "total"=total_ss, ss_results_summ)
    results_ss_summary_A<-rbind(results_ss_summary_A, t(week_results_ss))
    
    print("Results saved")
  }
}

saveRDS(results_A, "data/output_data/exp_A_bulk_analysis.RData")
saveRDS(results_ss_summary_A, "data/output_data/exp_A_bulk_analysis_subsampled.RData")

# repeat for beta

rm(exp_A_long)
rm(results_A)
rm(results_ss_summary_A)

min_depth<-min(exp_B_long$total_exp) # total contains the total number of sequences from each week
subsampleSize_b<-round(min_depth/10, 0)*10

results_B<-data.frame()
results_ss_summary_B<-data.frame()

for (id in unique(exp_B_long$ID)){
  print(paste("ID:", id))
  thisID<-exp_B_long[exp_B_long$ID == id,]
  sum(thisID$counts)
  for (week in unique(thisID$week_PCR)){
    print(paste("week:", week))
    thisWeek<-thisID[(thisID$week_PCR == week) & (thisID$counts > 0),]
    print(dim(thisWeek))
    control<-unique(thisWeek$control)
    total<-unique(thisWeek$total_exp)
    
    expanded<-thisWeek %>% uncount(counts)
    assert("Assert expected number of rows in table", total == dim(expanded)[1])
    
    myresults<-extract_diversity_metrics(thisWeek, expanded)
    
    week_results<-c("id"=id, "control" = control, "week"=week, "total"=total, 
                    myresults)
    results_B<-rbind(results_B, t(week_results))
    
    ss_results_all<-data.frame()
    
    for (i in 1:100){
      subsampled_expanded<-data.frame(expanded[sample(nrow(expanded), subsampleSize_b), ])
      total_ss<-subsampleSize_b
      
      subsampled<- subsampled_expanded %>% group_by(across(everything())) %>% count()
      subsampled$counts<-subsampled$n
      assert("Collapsing table correctly", dim(subsampled)[1] == dim(subsampled_expanded[!duplicated(subsampled_expanded),])[1])
      
      myresults_ss<-extract_diversity_metrics(subsampled, subsampled_expanded)
      ss_results_all<-rbind(ss_results_all, t(myresults_ss))
      
    }
    
    ss_results_summ<-colMeans(ss_results_all)
    
    week_results_ss<-c("id"=id, "control" = control, "week"=week, "total"=total_ss, ss_results_summ)
    results_ss_summary_B<-rbind(results_ss_summary_B, t(week_results_ss))
  }
}

saveRDS(results_B, "data/output_data/exp_B_bulk_analysis.RData")
saveRDS(results_ss_summary_B, "data/output_data/exp_B_bulk_analysis_subsampled.RData")
