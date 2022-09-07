#new script which makes a control set of TCRs equal to the expanded set for alpha and betza
#on Linux
drive <- "/media/benny/data/"
#on Windows
drive<-"C:/users/Benny Chain/"


#this is the path to the master files on my computer
input_data<-"Dropbox/Temp (1)/COVID-19/Data/"

file<-"exp_AB_wide3.RData"

file_name<-load(paste0(drive,input_data,file))
chain<-"alpha"
chain<-"beta"

i_chain<-which(exp_AB_wide1$chain==chain)
i_ctrl<-which(exp_AB_wide1$control==FALSE)
inter<-intersect(i_chain,i_ctrl)

#load controls

if (chain=="alpha"){
  load(file=paste0(drive,input_data,"all_A_long.RData"))
  data<-all_A_long
rm(all_A_long)}

if (chain=="beta"){
  load(file=paste0(drive,input_data,"all_B_long.RData"))
  data<-all_B_long
rm(all_B_long)}

#make a control set of the sam esize as expanded non-controls
#control_a<-list()
control_b<-list()

i<-1
for (i in 1:10){
  i_ss<-sample(1:length(unique(data$junction_aa)),length(inter))
  if (chain == "alpha"){control_a[[i]]<-unique(data$junction_aa)[i_ss]}
  if (chain == "beta"){control_b[[i]]<-unique(data$junction_aa)[i_ss]}
  }
  
  #save(control_a,file=paste0(drive,input_data,"control_alpha_1.RData"))
  
  save(control_b,file=paste0(drive,input_data,"control_beta_1.RData"))
  


