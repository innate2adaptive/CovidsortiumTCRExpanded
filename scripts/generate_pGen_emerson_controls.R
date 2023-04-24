# subset all_A_long and all_B_long to the controls I use 
# for the emerson analysis and for pGen

library(dplyr)

options(timeout = max(10000, getOption("timeout"))) # increasing timeout might be necessary to load the files

myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

# the data is in long format, so first I keep only one entry for each tcr in each patient (here I have multiple entries for multiple timepoints)
all_B_long1<-all_B_long[!(duplicated(all_B_long[, c("decombinator_id", "ID")])),]
rm(all_B_long)
# remove PCR- controls
all_B_long2<-all_B_long1[all_B_long1$control == FALSE,]
rm(all_B_long1)

load("data/output_data/exp_AB_wide4.RData")
exp_AB_wide3<-exp_AB_wide3[exp_AB_wide3$control == FALSE,]
exp_a<-exp_AB_wide3[exp_AB_wide3$chain == "alpha",]
exp_b<-exp_AB_wide3[exp_AB_wide3$chain == "beta",]

# remove expanded from controls

all_B_long2$ID<-as.character(all_B_long2$ID)
exp_b$ID<-as.character(exp_b$ID)
nonexp_B_long<-anti_join(all_B_long2, exp_b, by=c("decombinator_id", "ID")) # remove expanded from the list
rm(all_B_long2)
nonexp_B_long<-nonexp_B_long[!duplicated(nonexp_B_long[,c("v_call", "junction_aa")]),]
write.csv(nonexp_B_long, file=gzfile("data/output_data/HCW_B_controls_emerson_pGen.csv.gz"))
rm(nonexp_B_long)

# repeat for alpha

myURL<-"https://www.dropbox.com/s/a7ymcecnpomge2e/all_A_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

# the data is in long format, so first I keep only one entry for each tcr in each patient (here I have multiple entries for multiple timepoints)
all_A_long1<-all_A_long[!(duplicated(all_A_long[, c("decombinator_id", "ID")])),]
rm(all_A_long)
# remove PCR- controls
all_A_long2<-all_A_long1[all_A_long1$control == FALSE,]
rm(all_A_long1)

# remove expanded from controls

all_A_long2$ID<-as.character(all_A_long2$ID)
exp_a$ID<-as.character(exp_a$ID)
nonexp_A_long<-anti_join(all_A_long2, exp_a, by=c("decombinator_id", "ID")) # remove expanded from the list
rm(all_A_long2)
nonexp_A_long<-nonexp_A_long[!duplicated(nonexp_A_long[,c("v_call", "junction_aa")]),]
write.csv(nonexp_A_long, file=gzfile("data/output_data/HCW_A_controls_emerson_pGen.csv.gz"))
rm(nonexp_A_long)
