#reanalyse TCRs and remove all those very high in weeks -2 or -3

file2<-"data/output_data/exp_AB_wide.RData"

load(file2)
print(dim(exp_AB_wide))

#remove TCRs which are up early
i_up1<-which(exp_AB_wide["proportion_-3"]>4)
i_up2<-which(exp_AB_wide["proportion_-2"]>4)

i_up<-unique(c(i_up1,i_up2))
exp_AB_wide1<-exp_AB_wide[-i_up,]
print(dim(exp_AB_wide1))

save(exp_AB_wide1, file="data/output_data/exp_AB_wide1.RData")
