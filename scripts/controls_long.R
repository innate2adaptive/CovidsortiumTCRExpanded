#new script which makes a control set of TCRs equal to the expanded set for alpha and beta

file_name<-load("data/output_data/exp_AB_wide3.RData")

control_a<-list()
control_b<-list()

for (chain in c("alpha", "beta")){

i_chain<-which(exp_AB_wide3$chain==chain)
i_ctrl<-which(exp_AB_wide3$control==FALSE)
inter<-intersect(i_chain,i_ctrl)

options(timeout = max(1000, getOption("timeout"))) # increasing timeout might be necessary to load the files

#load controls

if (chain=="alpha"){
  myURL<-"https://www.dropbox.com/s/a7ymcecnpomge2e/all_A_long.RData?raw=1"
  myConnection <- url(myURL)
  print(load(myConnection))
  close(myConnection)
  data<-all_A_long
rm(all_A_long)}

if (chain=="beta"){
  myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
  myConnection <- url(myURL)
  print(load(myConnection))
  close(myConnection)
  data<-all_B_long
rm(all_B_long)}

#make a control set of the same size as expanded non-controls

i<-1
for (i in 1:10){
  i_ss<-sample(1:length(unique(data$junction_aa)),length(inter))
  if (chain == "alpha"){control_a[[i]]<-unique(data$junction_aa)[i_ss]}
  if (chain == "beta"){control_b[[i]]<-unique(data$junction_aa)[i_ss]}
  }
rm(data)

} # close loop over chains

save(control_a,file="data/output_data/control_alpha_1.RData")
save(control_b,file="data/output_data/control_beta_1.RData")
  


