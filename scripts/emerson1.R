# load and compute emerson sharing

sharing<-readRDS("data/downloaded_data/Emerson_sharing_aa.rds")
load("data/output_data/exp_AB_wide3.RData")

sharing_in_exp<-sharing[sharing$aminoAcid %in% exp_AB_wide3$junction_aa,]

myURL<-"https://www.dropbox.com/s/a7ymcecnpomge2e/all_A_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

non_exp_A<-all_A_long[!(all_A_long$junction_aa %in% unique(exp_AB_wide3$junction_aa)),]
sharing_in_ctrl_A<-sharing[sharing$aminoAcid %in% unique(non_exp_A$junction_aa),]
rm(all_A_long)
rm(non_exp_A)

myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

non_exp_B<-all_B_long[!(all_B_long$junction_aa %in% unique(exp_AB_wide3$junction_aa)),]
sharing_in_ctrl_B<-sharing[sharing$aminoAcid %in% unique(non_exp_B$junction_aa),]

rm(all_B_long)
rm(non_exp_B)
rm(exp_AB_wide3)
rm(sharing)

pX_0 <- exp(-sharing_in_exp$sharing_level/(10^5))
f<--log10(pX_0)*10

hist(sharing_in_exp$sharing_level)
hist(sharing_in_ctrl_A$sharing_level)
