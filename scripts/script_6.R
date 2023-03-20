"
This script takes the expanded and the non-expanded control sequences and plots the differences
in clustering across a number of different thresholds.
This is a little slow, be patient!
"
require(kernlab)
require(igraph)
require(ggplot2)
require(reshape2)

sk <- stringdot(type="spectrum", length=3, normalized=TRUE)

#folder for plots
folder_plots<-"output_figures/" # don't save to dropbox, so if people use this will save locally
#folder for data output
folder_data<-"data/output_data/"  # don't save to dropbox, so if people use this will save locally

# load data - using low p(gen) controls
load("data/output_data/exp_AB_wide3.RData")
load("data/output_data/control_alpha_1.RData") # this is a set of non-expanded sequences selected at random
load("data/output_data/control_beta_1.RData") # this is a set of non-expanded sequences selected at random 

exp_a<-exp_AB_wide3[exp_AB_wide3$chain=="alpha",]
exp_b<-exp_AB_wide3[exp_AB_wide3$chain=="beta",]
exp_a.noC<-exp_a[exp_a$control==FALSE,] # remove PCR-
exp_b.noC<-exp_b[exp_b$control==FALSE,] # remove PCR-

a_exp.noC<-unique(exp_a.noC$junction_aa)
b_exp.noC<-unique(exp_b.noC$junction_aa)

kma_exp.noC<-kernelMatrix(sk, a_exp.noC, a_exp.noC)
kmb_exp.noC<-kernelMatrix(sk, b_exp.noC, b_exp.noC)

### create km for non-expanded controls (x10 for each)

c.a<-list()
c.b<-list()

for (i in 1:10){
  print(i)
  kma.c<-kernelMatrix(sk, unique(unlist(control_a[i])), unique(unlist(control_a[i])))
  kmb.c<-kernelMatrix(sk, unique(unlist(control_b[i])), unique(unlist(control_b[i])))
  c.a[[i]]<-kma.c
  c.b[[i]]<-kmb.c
}

### Look at how the clustering performs at different thresholds
### various metrics for clustering at different thresholds

percentage_in_cluster <- function(clus.g){
  clus.g.count<-table(clus.g)
  clus.g.nosinglet<-clus.g.count[clus.g.count > 1]
  tot<-sum(clus.g.count)
  in.clus<-sum(clus.g.nosinglet)
  perc<-in.clus/tot
  return(perc)
}

percentage_in_large_cluster <- function(clus.g){
  clus.g.count<-table(clus.g)
  clus.g.nosinglet<-clus.g.count[clus.g.count > 4]
  tot<-sum(clus.g.count)
  in.clus<-sum(clus.g.nosinglet)
  perc<-in.clus/tot
  return(perc)
}

avg_clus_size<-function(clus.g){
  clus.g.count<-table(clus.g)
  clus.g.nosinglet<-clus.g.count[clus.g.count > 1]
  avg<-unlist(mean(clus.g.nosinglet))
  return(avg)
}

num_large_clusters <- function(clus.g){
  clus.g.count<-table(clus.g)
  large.clus<-length(clus.g.count[clus.g.count > 4])
  return(large.clus)
}


tha<-seq(0, 1, by=0.02)
thb<-seq(0, 1, by=0.02)

clus.size.count.a<-list()
clus.size.count.a.noC<-list()
clus.size.count.a.c<-list()
clus.size.count.b<-list()
clus.size.count.b.noC<-list()
clus.size.count.b.c<-list()

perc.clustering.a<-list()
perc.clustering.large.a<-list()
avg.cluster.size.a<-list()
large.clusters.a<-list()

perc.clustering.b<-list()
perc.clustering.large.b<-list()
avg.cluster.size.b<-list()
large.clusters.b<-list()

perc.clustering.a.noC<-list()
perc.clustering.large.a.noC<-list()
avg.cluster.size.a.noC<-list()
large.clusters.a.noC<-list()

perc.clustering.b.noC<-list()
perc.clustering.large.b.noC<-list()
avg.cluster.size.b.noC<-list()
large.clusters.b.noC<-list()

perc.clustering.a.c<-list()
perc.clustering.large.a.c<-list()
avg.cluster.size.a.c<-list()
large.clusters.a.c<-list()

perc.clustering.b.c<-list()
perc.clustering.large.b.c<-list()
avg.cluster.size.b.c<-list()
large.clusters.b.c<-list()

for (j in 1:length(tha)){ # iterate over various thresholds
  ta<-tha[j]
  tb<-thb[j]
  print(paste("threshold:", ta))
  
  clus.size.count.a.c[[j]]<-list()
  perc.clustering.a.c[[j]]<-list()
  perc.clustering.large.a.c[[j]]<-list()
  avg.cluster.size.a.c[[j]]<-list()
  large.clusters.a.c[[j]]<-list()

  clus.size.count.b.c[[j]]<-list()
  perc.clustering.b.c[[j]]<-list()
  perc.clustering.large.b.c[[j]]<-list()
  avg.cluster.size.b.c[[j]]<-list()
  large.clusters.b.c[[j]]<-list()

  print("2")
#
  connect_a.noC<-(kma_exp.noC>=ta)
  km_adj.a.noC<-kma_exp.noC
  km_adj.a.noC[which(connect_a.noC)]<-1
  km_adj.a.noC[which(!connect_a.noC)]<-0

  connect_b.noC<-(kmb_exp.noC>=tb)
  km_adj.b.noC<-kmb_exp.noC
  km_adj.b.noC[which(connect_b.noC)]<-1
  km_adj.b.noC[which(!connect_b.noC)]<-0

  print("5")
#
  km_graph.a.noC<-graph.adjacency(km_adj.a.noC, mode="undirected",weighted = TRUE,add.colnames=NULL,diag=FALSE)
  clus.a.noC<-clusters(km_graph.a.noC)$membership
  clus.a.noC.size.count<-data.frame(table(table(clus.a.noC)))
  perc.a.noC<-percentage_in_cluster(clus.a.noC)
  perc.large.a.noC<-percentage_in_large_cluster(clus.a.noC)
  avg.clus.a.noC<-avg_clus_size(clus.a.noC)
  num.large.clus.a.noC<-num_large_clusters(clus.a.noC)
  clus.size.count.a.noC<-c(clus.size.count.a.noC, clus.a.noC.size.count)
  perc.clustering.a.noC<-c(perc.clustering.a.noC, perc.a.noC)
  perc.clustering.large.a.noC<-c(perc.clustering.large.a.noC, perc.large.a.noC)
  avg.cluster.size.a.noC<-c(avg.cluster.size.a.noC, avg.clus.a.noC)
  large.clusters.a.noC<-c(large.clusters.a.noC, num.large.clus.a.noC)

  print("6")

  km_graph.b.noC<-graph.adjacency(km_adj.b.noC, mode="undirected",weighted = TRUE,add.colnames=NULL,diag=FALSE)
  clus.b.noC<-clusters(km_graph.b.noC)$membership
  clus.b.noC.size.count<-data.frame(table(table(clus.b.noC)))
  perc.b.noC<-percentage_in_cluster(clus.b.noC)
  perc.large.b.noC<-percentage_in_large_cluster(clus.b.noC)
  avg.clus.b.noC<-avg_clus_size(clus.b.noC)
  num.large.clus.b.noC<-num_large_clusters(clus.b.noC)
  clus.size.count.b.noC<-c(clus.size.count.b.noC, clus.b.noC.size.count)
  perc.clustering.b.noC<-c(perc.clustering.b.noC, perc.b.noC)
  perc.clustering.large.b.noC<-c(perc.clustering.large.b.noC, perc.large.b.noC)
  avg.cluster.size.b.noC<-c(avg.cluster.size.b.noC, avg.clus.b.noC)
  large.clusters.b.noC<-c(large.clusters.b.noC, num.large.clus.b.noC)

  print("now controls")

  for (i in 1:10){ # 10 controls for each thresh
    print(paste0("controls, ", i))
    kma.c<-c.a[[i]]
    connect_a<-(kma.c>=ta)
    km_adj.a.c<-kma.c
    km_adj.a.c[which(connect_a)]<-1
    km_adj.a.c[which(!connect_a)]<-0
  #
    print("c2")
  #
    kmb.c<-c.b[[i]]
    connect_b<-(kmb.c>=tb)
    km_adj.b.c<-kmb.c
    km_adj.b.c[which(connect_b)]<-1
    km_adj.b.c[which(!connect_b)]<-0
  #
    print("c3")
  #
    km_graph.a.c<-graph.adjacency(km_adj.a.c, mode="undirected",weighted = TRUE,add.colnames=NULL,diag=FALSE)
    clus.a.c<-clusters(km_graph.a.c)$membership
    clus.a.c.size.count<-data.frame(table(table(clus.a.c)))
    perc.a.c<-percentage_in_cluster(clus.a.c)
    perc.large.a.c<-percentage_in_large_cluster(clus.a.c)
    avg.clus.a.c<-avg_clus_size(clus.a.c)
    num.large.clus.a.c<-num_large_clusters(clus.a.c)
    clus.size.count.a.c[[j]]<-c(clus.size.count.a.c[[j]], clus.a.c.size.count)
    perc.clustering.a.c[[j]]<-c(perc.clustering.a.c[[j]], perc.a.c)
    perc.clustering.large.a.c[[j]]<-c(perc.clustering.large.a.c[[j]], perc.large.a.c)
    avg.cluster.size.a.c[[j]]<-c(avg.cluster.size.a.c[[j]], avg.clus.a.c)
    large.clusters.a.c[[j]]<-c(large.clusters.a.c[[j]], num.large.clus.a.c)
  #
    print("c4")
  #
    km_graph.b.c<-graph.adjacency(km_adj.b.c, mode="undirected",weighted = TRUE,add.colnames=NULL,diag=FALSE)
    clus.b.c<-clusters(km_graph.b.c)$membership
    clus.b.c.size.count<-data.frame(table(table(clus.b.c)))
    perc.b.c<-percentage_in_cluster(clus.b.c)
    perc.large.b.c<-percentage_in_large_cluster(clus.b.c)
    avg.clus.b.c<-avg_clus_size(clus.b.c)
    num.large.clus.b.c<-num_large_clusters(clus.b.c)
    clus.size.count.b.c[[j]]<-c(clus.size.count.b.c[[j]], clus.b.c.size.count)
    perc.clustering.b.c[[j]]<-c(perc.clustering.b.c[[j]], perc.b.c)
    perc.clustering.large.b.c[[j]]<-c(perc.clustering.large.b.c[[j]], perc.large.b.c)
    avg.cluster.size.b.c[[j]]<-c(avg.cluster.size.b.c[[j]], avg.clus.b.c)
    large.clusters.b.c[[j]]<-c(large.clusters.b.c[[j]], num.large.clus.b.c)
  }
}

### functions to rearrange results in more useful format
rearrange_clus_size<-function(x, th){
  control_counts<-data.frame()
  for (j in seq(2, 2*length(th), 2)){
    count<-x[[j]]
    # print(count)
    c.size<-as.numeric(levels(x[[j-1]]))
    # print(c.size)
    xdf<-data.frame(cbind(count, c.size))
    xdf["type"]<-"data"
    xdf["thresh"]<-th[j/2]
    # print(xdf)
    control_counts<-rbind(control_counts, xdf)
  }
  return(control_counts)
}
rearrange_clus_size_ctrls<-function(x, th){
  control_counts<-data.frame()
  for (j in seq(1, length(th))){
    for (i in seq(2,20,2)){
      count<-x[[j]][[i]]
      c.size<-as.numeric(levels(x[[j]][[i-1]]))
      xdf<-data.frame(cbind(count, c.size))
      xdf["type"]<-paste("control", i/2, sep="")
      xdf["thresh"]<-th[j]
      control_counts<-rbind(control_counts, xdf)
    }
  }
  return(control_counts)
}

## clean data

large.clus.a.c<-data.frame(do.call(rbind, large.clusters.a.c))
large.clus.a<-data.frame(do.call(rbind, large.clusters.a.noC))
large.clus.a$thres<-tha
large.clus.a.c$thres<-tha
avg.clus.a.c<-data.frame(do.call(rbind,avg.cluster.size.a.c))
avg.clus.a<-data.frame(do.call(rbind,avg.cluster.size.a.noC))
avg.clus.a.c$thres<-tha
avg.clus.a$thres<-tha
perc.clus.a.c<-data.frame(do.call(rbind,perc.clustering.a.c))
perc.clus.a<-data.frame(do.call(rbind,perc.clustering.a.noC))
perc.clus.a.c$thres<-tha
perc.clus.a$thres<-tha
perc.clus.large.a.c<-data.frame(do.call(rbind,perc.clustering.large.a.c))
perc.clus.large.a<-data.frame(do.call(rbind,perc.clustering.large.a.noC))
perc.clus.large.a.c$thres<-tha
perc.clus.large.a$thres<-tha
clus.size.count.a.c.df<-rearrange_clus_size_ctrls(clus.size.count.a.c,tha)
clus.size.count.a.df<-rearrange_clus_size(clus.size.count.a.noC,tha)
clus.size.count.a.all<-rbind(clus.size.count.a.c.df, clus.size.count.a.df)
clus.size.count.a.all$tot.count<-clus.size.count.a.all$c.size*clus.size.count.a.all$count
clus.size.count.a.all$prop<-ifelse(clus.size.count.a.all$type=="data", clus.size.count.a.all$tot.count/dim(exp_a)[1], clus.size.count.a.all$tot.count/length(control_a[[1]]))
clus.size.count.a.all$c.size.cat<-ifelse(clus.size.count.a.all$c.size == 1, "singlet",
                                         ifelse(clus.size.count.a.all$c.size == 2, "doublet",
                                                ifelse(clus.size.count.a.all$c.size < 5, "3-4",
                                                       ifelse(clus.size.count.a.all$c.size < 11, "5-10",
                                                              ifelse(clus.size.count.a.all$c.size < 101, "11-100", ">100")))))
mylevels<-c(">100", "11-100", "5-10", "3-4", "doublet", "singlet")
clus.size.count.a.all$c.size.cat<-factor(clus.size.count.a.all$c.size.cat, levels=mylevels)
cols<-c("chocolate4", "lightsalmon", "moccasin", "lightskyblue", "lightsteelblue1", "lightsteelblue")
names(cols)<-mylevels
 
large.clus.b.c<-data.frame(do.call(rbind, large.clusters.b.c))
large.clus.b<-data.frame(do.call(rbind, large.clusters.b.noC))
large.clus.b$thres<-thb
large.clus.b.c$thres<-thb
avg.clus.b.c<-data.frame(do.call(rbind,avg.cluster.size.b.c))
avg.clus.b<-data.frame(do.call(rbind,avg.cluster.size.b.noC))
avg.clus.b.c$thres<-thb
avg.clus.b$thres<-thb
perc.clus.b.c<-data.frame(do.call(rbind,perc.clustering.b.c))
perc.clus.b<-data.frame(do.call(rbind,perc.clustering.b.noC))
perc.clus.b.c$thres<-thb
perc.clus.b$thres<-thb
perc.clus.large.b.c<-data.frame(do.call(rbind,perc.clustering.large.b.c))
perc.clus.large.b<-data.frame(do.call(rbind,perc.clustering.large.b.noC))
perc.clus.large.b.c$thres<-thb
perc.clus.large.b$thres<-thb
clus.size.count.b.c.df<-rearrange_clus_size_ctrls(clus.size.count.b.c,thb)
clus.size.count.b.df<-rearrange_clus_size(clus.size.count.b.noC,thb)
clus.size.count.b.all<-rbind(clus.size.count.b.c.df, clus.size.count.b.df)
clus.size.count.b.all$tot.count<-clus.size.count.b.all$c.size*clus.size.count.b.all$count
clus.size.count.b.all$prop<-ifelse(clus.size.count.b.all$type=="data", clus.size.count.b.all$tot.count/dim(exp_b)[1], clus.size.count.b.all$tot.count/length(control_b[[1]]))
clus.size.count.b.all$c.size.cat<-ifelse(clus.size.count.b.all$c.size == 1, "singlet",
                                         ifelse(clus.size.count.b.all$c.size == 2, "doublet",
                                                ifelse(clus.size.count.b.all$c.size < 5, "3-4",
                                                       ifelse(clus.size.count.b.all$c.size < 11, "5-10",
                                                              ifelse(clus.size.count.b.all$c.size < 101, "11-100", ">100")))))
clus.size.count.b.all$c.size.cat<-factor(clus.size.count.b.all$c.size.cat, levels=mylevels)

## more functions to faciliate plotting

prepare_results<-function(x, x1){
  x<-data.frame(x)
  for (i in 1:11){
    y<-names(x)[i]
    x[,y]<-as.numeric(x[,y])
  }
  # print(sapply(x, class))
  avgs<-sapply(data.frame(t(x[,!names(x)=="thres"])), mean)
  stds<-sapply(data.frame(t(x[,!names(x)=="thres"])), sd)
  x$avg<-as.numeric(avgs)
  x$std<-as.numeric(stds)
  x.m<-melt(x, id.vars=c("thres", "avg", "std"))

  x1<-data.frame(x1)
  names(x1)[1]<-"data"
  for (i in 1:2){
    y<-names(x1)[i]
    x1[,y]<-as.numeric(x1[,y])
  }
  x1$avg<-NA
  x1$std<-NA

  x1.m<-melt(x1, id.vars=c("thres", "avg", "std"))
  x0<-rbind(x.m, x1.m)

  return(x0)
}

remove_dup_na<-function(x){
  y<-x[!duplicated(x[c("thres", "avg")]),]
  x1<-y[!is.na(y[c("avg")]),]
  return(x1)
}

plot_results<-function(x, shape=16){
  mycolors<-c("grey", "grey","grey","grey","grey","grey",
              "grey","grey","grey","grey", "red")
  p<- ggplot(x, aes(x=as.numeric(thres), y=as.numeric(value), color=variable)) +
    geom_line() + geom_point(size=2, shape=shape) +
    geom_line(aes(x=thres, y=as.numeric(avg)), color="black",
              data=function(x){remove_dup_na(x)}) +
    geom_point(aes(x=thres, y=as.numeric(avg)), color="black",
               data=function(x){remove_dup_na(x)}, size=2, shape=shape) +
    geom_errorbar(aes(ymin=avg-std, ymax=avg+std), color="black", width=0.005,
                  data=function(x){remove_dup_na(x)}) +
    scale_color_manual(values = mycolors, guide = "none")
    # ylim(0,1)
  return(p)
}

## plot all

perc.a<-prepare_results(perc.clus.a.c, perc.clus.a)
p.perc.a<-plot_results(perc.a)
pa <-p.perc.a +
  geom_vline(xintercept = 0.76, lty="dashed") +
  xlab("threshold set for clustering") +
  ylab("Proportion clustering") +
  ggtitle("Proportion of CDR3s in clusters, alpha") + 
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=16),
        title=element_text(size=14))

svg("output_figures/FigS11_alpha.svg")
print(pa)
dev.off()

perc.b<-prepare_results(perc.clus.b.c, perc.clus.b)
p.perc.b<-plot_results(perc.b)
pb<-p.perc.b +
  geom_vline(xintercept = 0.72, lty="dashed") +
  xlab("threshold set for clustering") +
  ylab("Proportion clustering") +
  ggtitle("Proportion of CDR3s in clusters, beta") +
  theme_classic() + # keep the theme consistent for all our plots
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=16),
        title=element_text(size=14))

svg("output_figures/FigS11_beta.svg")
print(pb)
dev.off()
