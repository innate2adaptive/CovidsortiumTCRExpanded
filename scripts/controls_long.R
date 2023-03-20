#new script which makes a control set of TCRs equal to the expanded set for alpha and beta
library(dplyr)

file_name<-load("data/output_data/exp_AB_wide3.RData")

control_a<-list()
control_b<-list()

for (chain in c("alpha", "beta")){

i_chain<-which(exp_AB_wide3$chain==chain)
i_ctrl<-which(exp_AB_wide3$control==FALSE)
inter<-intersect(i_chain,i_ctrl)
exp<-exp_AB_wide3[inter,]

options(timeout = max(1000, getOption("timeout"))) # increasing timeout might be necessary to load the files

#load controls

if (chain=="alpha"){
  myURL<-"https://www.dropbox.com/s/a7ymcecnpomge2e/all_A_long.RData?raw=1"
  myConnection <- url(myURL)
  print(load(myConnection))
  close(myConnection)
  all_A_long1<-all_A_long[!(duplicated(all_A_long[, c("decombinator_id", "ID")])),]
  all_A_long1<-all_A_long1[all_A_long1$control == FALSE,]
  data0<-all_A_long1
  rm(all_A_long)
  rm(all_A_long1)
  #remove MAIT
  i_V<-which(data0$v_call=="TRAV1-2")
  i_J<-which(data0$j_call=="TRAJ12" | data0$j_call=="TRAJ20" | data0$j_call=="TRAJ33")
  data1<-data0[-intersect(i_V,i_J),]
  
  #remove iKT
  i_V<-which(data1$v_call=="TRAV10")
  i_J<-which(data1$j_call=="TRAJ18")
  data2<-data1[-(intersect(i_V,i_J)),]

  nonexp_A_long<-anti_join(data2, exp, by=c("decombinator_id", "ID")) # takes non-expanded from PCR+
  data2<-nonexp_A_long

rm(data0)
rm(data1)}

if (chain=="beta"){
  myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
  myConnection <- url(myURL)
  print(load(myConnection))
  close(myConnection)
  all_B_long1<-all_B_long[!(duplicated(all_B_long[, c("decombinator_id", "ID")])),]
  all_B_long1<-all_B_long1[all_B_long1$control == FALSE,]
  nonexp_B_long<-anti_join(all_B_long1, exp, by=c("decombinator_id", "ID")) # takes non-expanded from PCR+
  data2<-nonexp_B_long
rm(all_B_long)
rm(all_B_long1)}

#make a control set of the same size as expanded non-controls

i<-1
for (i in 1:10){
  i_ss<-sample(1:length(unique(data2$junction_aa)),length(inter))
  if (chain == "alpha"){control_a[[i]]<-unique(data2$junction_aa)[i_ss]}
  if (chain == "beta"){control_b[[i]]<-unique(data2$junction_aa)[i_ss]}
  }
rm(data)

} # close loop over chains

save(control_a,file="data/output_data/control_alpha_1.RData")
save(control_b,file="data/output_data/control_beta_1.RData")
  


