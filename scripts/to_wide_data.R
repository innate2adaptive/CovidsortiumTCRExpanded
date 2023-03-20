#look at data from COVID19 using Marc's master files
#plot COVID response
library(mgcv)
library(tidyr)
library(dplyr)
library(ggplot2)

#read dates of PCR positivity
PCR_pos2<-read.table(file="data/PCR_positive.csv",header=TRUE,stringsAsFactors = FALSE,colClasses=c("character", "numeric"),sep=",")
#remove 241 which does not exist in TCR data
PCR_pos1<-PCR_pos2[-27,]
PCR_pos<-PCR_pos1$Week.first..
names(PCR_pos)<-as.character(PCR_pos1[,1])
rm(PCR_pos2)
#PCR_pos

options(timeout = max(1000, getOption("timeout")))

for (chain in c("alpha", "beta")){
#load total and expanded data either for alpha and beta. Combine at the end.
if (chain=="alpha"){
#load all data master file
  myURL <- "https://www.dropbox.com/s/93jzrddo291h202/data_tsv_alpha.Rdata?raw=1"
  myConnection <- url(myURL)
  print(load(myConnection))
  close(myConnection)
#load expanded TCRs
  load( file = "data/output_data/TCR_change_all_alpha.RData")
  TCR_change_HCW<-TCR_change_HCW_a
  data<-data_tsv_alpha
  rm(data_tsv_alpha)
  rm(TCR_change_HCW_a)
}

#Beta
if (chain=="beta"){
  myURL<-"https://www.dropbox.com/s/7osyel4pit0gayn/data_tsv_beta.Rdata?raw=1"
  myConnection <- url(myURL)
  print(load(myConnection))
  close(myConnection)
#load data on expanded TCRs
load( file = "data/output_data/TCR_change_all_beta.RData")
  TCR_change_HCW<-TCR_change_HCW_b
  data<-data_tsv_beta
  rm(data_tsv_beta)
  rm(TCR_change_HCW_b)
}

#convert time points to numeric and add a column to data
tp_week<-rep(100,dim(data)[1])
tp_week[which(data$timepoint=="BL")]<-0
tp_split<-strsplit(data$timepoint,"U")
i_spl<-which(sapply(tp_split,length)==2)
week<-as.numeric(sapply(tp_split[i_spl],function(x){x[2]}))
tp_week[i_spl]<-week
data$tp_week<-tp_week

########################################################################################
####################################################################################################
#cycle through and extract all the changed TCRs for all time points from the big data file
i<-2

#this file collects all the expanded TCRs
data_exp_all<-c()
  for (i in 1 : length(TCR_change_HCW)){
  #for (i in 1 : 4){
  ID<-names(TCR_change_HCW)[[i]]
  if(!is.na(TCR_change_HCW[[i]][1])){
  TCR_exp<-unique(rownames(TCR_change_HCW[[i]]))
  #select from master file just TCRs from this individual
  data_id<-data[which(data$ID==ID),]
  #add a column with total number of TCRs so can calculate proportions later
  #data_id$total<-sum(data_id$duplicate_count)
  data_id$total<- as.integer(round(data_id$duplicate_count/data_id$proportion,0))
  #remove the character descprtions of time because redundant; also remove proportion
  data_id<-subset(data_id, select=-c(timepoint,proportion))

  #combine weeks and totals into one field so can pivot
  data_id$week_total<-paste(data_id$tp_week,data_id$total,sep="_")
  #add a column with week which became PCR +
  data_id$PCR<-PCR_pos[as.character(as.numeric(as.character(ID)))]
  #data_id$PCR[1]
  #if(length(i_ID)>0){data_id$PCR<-PCR_pos1[i_ID,2]} else {data_id$PCR<-NA}
  #available times and totals for this ID
  tp_week1<-unique(data_id$week_total)
  #timepoints<-unique(data_id$timepoint)
  #names_time<-paste(timepoints,tp_week,sep="&")
  #construct a unique TCR identifier
  data_id_TCR<-paste(data_id$junction,data_id$v_call,data_id$j_call,data_id$decombinator_id,  sep="_" )

  #extract the changing TCRs from this
  TCR_exp_i<-which(data_id_TCR %in% TCR_exp)
  data_id_TCR_exp<-data_id[TCR_exp_i,]

  #a<-TCR_change_HCW[[1]]
  #go to wide form to incorporate all available timepoints
  #go to wide form to incorporate all available timepoints for this individual. Set values to 0 if absent.
  data_wide <- pivot_wider(data_id_TCR_exp,id_cols=c(1,2,3,4,5,c("ID","control","PCR")),names_from = c("week_total"),values_from = c("duplicate_count"), values_fill = 0)
  data_long<-pivot_longer(data_wide,cols = all_of(tp_week1), names_to = "tp_week", values_to ="duplicate_count")
  data_long$total<-matrix(as.numeric(unlist(strsplit(data_long$tp_week,split="_"))),ncol=2,byrow = TRUE)[,2]
  data_long$tp_week<-matrix(as.numeric(unlist(strsplit(data_long$tp_week,split="_"))),ncol=2,byrow = TRUE)[,1]
  data_long$chain<-chain

  ############################################################################################
  #make some further edits to this master file
  #add a column with proportion as TCRs per million
  TCR_million<-(10^6)*data_long$duplicate_count/data_long$total
  data_long$proportion<-TCR_million
  #renumber weeks so relative to covid PCR +
  if(data_long$control[1]=="FALSE"){week_PCR<-as.numeric(as.character(data_long$tp_week))-as.numeric(as.character(data_long$PCR))}
  if(data_long$control[1]=="TRUE"){week_PCR<-as.numeric(as.character(data_long$tp_week))}

  #rename the weeks so everything >= 12 is 14 and add this as an extra week
  week_PCR1<-week_PCR
  week_PCR1[week_PCR>11]<-14
  data_long$week_PCR<-week_PCR1
  
  print(ID)
  print(data_long[1:10,])

  #remove tp_week and duplicate_count columns because they duplicate the PCR_week and proportion columns
  data_long<-subset(data_long,select=-c(tp_week,duplicate_count))

data_exp_all<-rbind(data_exp_all,data_long)
cat(ID,"\n")
    }#end of if clause for those with only one time point
  }


#################################################################################
#rename the expanded TCRs in long and wide format, combine alpha and beta and save
if (chain=="alpha"){
  exp_A_long<-data_exp_all
  #save(exp_A_long,file=paste0(drive,input_data,"exp_A_long.RData"))
  #convert to wide format
  exp_A_wide <- pivot_wider(exp_A_long, id_cols=c(1,2,3,4,5,c("ID","control","PCR","chain")),names_from = c("week_PCR"),values_from = c("proportion","total"),names_sort = TRUE)
  #save this as the wide version of the expanded set of TCRS
  #save(exp_A_wide,file=paste0(drive,input_data,"exp_A_wide.RData"))
}

#rename the expanded TCRs in long and wide format
if (chain=="beta"){
  exp_B_long<-data_exp_all
  #save(exp_A_long,file=paste0(drive,input_data,"exp_A_long.RData"))
  #convert to wide format
  exp_B_wide <- pivot_wider(exp_B_long, id_cols=c(1,2,3,4,5,c("ID","control","PCR","chain")),names_from = c("week_PCR"),values_from = c("proportion","total"),names_sort = TRUE)
  #save this as the wide version of the expanded set of TCRS
  #save(exp_A_wide,file=paste0(drive,input_data,"exp_A_wide.RData"))
}
} # closes the loop over the two chains

exp_AB_long<-rbind(exp_A_long,exp_B_long)
exp_AB_wide<-rbind(exp_A_wide,exp_B_wide)

save(exp_AB_long,file="data/output_data/exp_AB_long.RData")
write.csv(exp_AB_long,file="data/output_data/exp_AB_long.csv")

save(exp_AB_wide,file="data/output_data/exp_AB_wide.RData")
write.csv(exp_AB_wide,file="data/output_data/exp_AB_wide.csv")
############################################################################
#######################################################################
