#Annotation of COVID-19 expanded TCRs
library(kernlab)
library(igraph)
library(stringdist)
library(tidyr)
library(dplyr)
library(data.table)
library(pheatmap)
library(ggplot2)
library(scales)

options(timeout = max(1000, getOption("timeout"))) # increasing timeout might be necessary to load the files
folder_plots<-"output_figures/"

myURL<-"https://www.dropbox.com/s/a7ymcecnpomge2e/all_A_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)
#plot number of reads per sample
abundance<-all_A_long$total*all_A_long$proportion/1E6
sample_ID<-paste(all_A_long$week_PCR,all_A_long$ID,sep="_")
sample_total<-aggregate(abundance,by=list(sample_ID),sum )

rm(all_A_long)
#betas
myURL<-"https://www.dropbox.com/s/9kl9s4y775wam9z/all_B_long.RData?raw=1"
myConnection <- url(myURL)
print(load(myConnection))
close(myConnection)
#plot number of reads per sample
abundance_B<-all_B_long$total*all_B_long$proportion/1E6
sample_ID_B<-paste(all_B_long$week_PCR,all_B_long$ID,sep="_")
sample_total_B<-aggregate(abundance_B,by=list(sample_ID_B),sum )
rm(all_B_long)

sample_total_A1<-sample_total
sample_total_A1$chain<-"alpha"

sample_total_B1<-sample_total_B
sample_total_B1$chain<-"beta"

all_counts<-rbind(sample_total_A1, sample_total_B1)
median(all_counts[all_counts$chain == "alpha",]$x)

median(all_counts[all_counts$chain == "beta",]$x)

p<-ggplot(all_counts) +
  geom_jitter(aes(x = chain, y = x), width = 0.1) +
  geom_boxplot(aes(x = chain, y = x), fill = NA, outlier.shape = NA) +
  labs(x = "", y = "read counts") +
  scale_y_log10(breaks = c(10^3, 10^4, 10^5, 10^6),
                labels = number_format(),
                limits = c(10^3, 10^6)) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        title=element_text(size=20)) 

svg(paste0(folder_plots, "FigS3A.svg"))
print(p)
dev.off()

rownames(sample_total)<-sample_total$Group.1
rownames(sample_total_B)<-sample_total_B$Group.1

paired_counts <- merge(sample_total, sample_total_B, by = "Group.1")

p1<-ggplot(paired_counts) +
  geom_point(aes(x = x.x, y = x.y), size = 2) +
  labs(x = "TCR alpha depth", y = "TCR beta depth") + 
  scale_y_log10(breaks = c(10^4, 10^4, 10^5),
                labels = number_format(),
                limits = c(10^4, 10^5.4)) +
  scale_x_log10(breaks = c(10^4, 10^4, 10^5),
                labels = number_format(),
                limits = c(10^4, 10^5.4)) +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        title=element_text(size=20)) 

svg(paste0(folder_plots, "FigS3B.svg"))
print(p1)
dev.off()

cor.test(paired_counts$x.x,paired_counts$x.y, method = "spearman")
