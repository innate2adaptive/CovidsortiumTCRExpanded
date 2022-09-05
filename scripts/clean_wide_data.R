#reanalyse TCRs and remove all those very high in weeks 2 or 3

#on Linux
drive <- "/media/benny/data/"
#on Windows
drive<-"C:/"


#this is the path to the master files on my computer
input_data<-"Dropbox/Temp (1)/COVID-19/Data/"

#working folder for collecting output
folder<-"Dropbox/R_temp/02_09_2021/"
dir.create(paste0(drive,folder))

file2<-"exp_AB_wide.RData"

load(paste0(drive,input_data,file1))
load(paste0(drive,input_data,file2))
#convert IDs to characters
exp_AB_wide_pos[,"ID"]<-exp_AB_wide$ID
#save again
save(exp_AB_wide_pos, file=paste0(drive,input_data,"exp_AB_wide_pos.RData"))
#remove TCRs which are up early
i_up1<-which(exp_AB_wide$`proportion_-3`>4)
i_up2<-which(exp_AB_wide$'proportion_-2'>4)

i_up<-unique(c(i_up1,i_up2))
exp_AB_wide1<-exp_AB_wide[-i_up,]

save(exp_AB_wide1, file=paste0(drive,input_data,"exp_AB_wide1.RData"))
