#Compare annotation of COVID-19 expanded and non-expanded TCRs
library(kernlab)
library(igraph)
library(stringdist)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)


#folder for plots
folder_plots<-"output_figures/" # don't save to dropbox, so if people use this will save locally
#folder for data output
folder_data<-"data/output_data/" 

file<-"data/output_data/exp_AB_wide3.RData"
file_name<-load(file)

options(timeout = max(1000, getOption("timeout")))

myURL<-"https://www.dropbox.com/s/a7ymcecnpomge2e/all_A_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

all_A_long1<-all_A_long[!(duplicated(all_A_long[, c("decombinator_id", "ID")])),]
rm(all_A_long)

myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)

all_B_long1<-all_B_long[!(duplicated(all_B_long[, c("decombinator_id", "ID")])),]
rm(all_B_long)

all_long <- rbind(all_A_long1, all_B_long1)
# save intermediate so I don't have to run this again and I can start over from here
save(all_long, "data/downloaded_data/all_long.RData") 

nonexp_long<-anti_join(all_long, exp_AB_wide3, by=c("decombinator_id", "ID"))