"
This script makes clustering plots for expanded and control sequences. It colors the clustering according to a number of metrics.
Note that the plots are saved out to an output folder, which is hard-coded.
"

require(kernlab)
require(igraph)
require(ggplot2)
require(reshape2)
library(qgraph)
library(RColorBrewer)
library(fields)
library(testit)
require(msa)
require(ggmsa)
# library(seqLogo)
require(ggseqlogo)

tha<-0.76
thb<-0.72
sk <- stringdot(type="spectrum", length=3, normalized=TRUE)

#folder for plots
folder_plots<-"output_figures/" # don't save to dropbox, so if people use this will save locally
#folder for data output
folder_data<-"data/output_data/"  # don't save to dropbox, so if people use this will save locally

load(paste0(folder_data, "exp_AB_wide3.RData")) # use the df without mait and invariant

exp.a<-exp_AB_wide3[(exp_AB_wide3$control==FALSE)&(exp_AB_wide3$chain=="alpha"),]
exp.b<-exp_AB_wide3[(exp_AB_wide3$control==FALSE)&(exp_AB_wide3$chain=="beta"),]

# exp.a.wc<-exp_AB_wide3[exp_AB_wide3$chain=="alpha",]
# exp.b.wc<-exp_AB_wide3[exp_AB_wide3$chain=="beta",]

# load annotated sequences
fs<-read.table("data/Francis_Science_CDR3.txt", header=TRUE)
tao<-read.table("data/Tao_CDR3.txt", header=TRUE)
load("data/VDJdb_single_specificity.RData")
vdj.cov<-VDJdb_all2[!is.na(VDJdb_all2$`SARS-CoV-2`),]

fs.a<-fs[fs$chain=="alpha",]
fs.b<-fs[fs$chain=="beta",]
tao.a<-tao[tao$gene=="alpha",]
tao.b<-tao[tao$gene=="beta",]
vdj.cov.a<-vdj.cov[vdj.cov$Gene=="TRA",]
vdj.cov.b<-vdj.cov[vdj.cov$Gene=="TRB",]

all.ann.a<-unique(c(fs.a$CDR3, tao.a$CDR3, vdj.cov.a$CDR3))
all.ann.b<-unique(c(fs.b$CDR3, tao.b$CDR3, vdj.cov.b$CDR3))

## load table where epitope names have been cleaned up
ann<-read.table("data/exp_AB_merge_2.txt", header=TRUE)

exp.a$ann<-unlist(lapply(exp.a$junction_aa, function(x) ifelse(x %in% all.ann.a, "ann", "not.ann")))
exp.b$ann<-unlist(lapply(exp.b$junction_aa, function(x) ifelse(x %in% all.ann.b, "ann", "not.ann")))

exp.a$cd4_cd8<-unlist(lapply(exp.a$junction_aa, function(x) ifelse(x %in% all.ann.a, toString(unique(ann[ann$junction_aa==x,]$CD8_CD4)), "na")))
exp.b$cd4_cd8<-unlist(lapply(exp.b$junction_aa, function(x) ifelse(x %in% all.ann.b, toString(unique(ann[ann$junction_aa==x,]$CD8_CD4)), "na")))

exp.a$public<-unlist(lapply(exp.a$junction_aa, function(x) length(unique(exp.a[exp.a$junction_aa == x, ]$ID))))
exp.b$public<-unlist(lapply(exp.b$junction_aa, function(x) length(unique(exp.b[exp.b$junction_aa == x, ]$ID))))

# exp.a.wc$ann<-unlist(lapply(exp.a.wc$junction_aa, function(x) ifelse(x %in% all.ann.a, "ann", "not.ann")))
# exp.b.wc$ann<-unlist(lapply(exp.b.wc$junction_aa, function(x) ifelse(x %in% all.ann.b, "ann", "not.ann")))

# assert("The CDR3s are the same", sort(unique(ann$CDR3)) == sort(unique(c(exp.a.wc[exp.a.wc$ann=="ann",]$junction_aa, exp.b.wc[exp.b.wc$ann=="ann",]$junction_aa))))

exp.a$pep<-unlist(lapply(exp.a$junction_aa, function(x) {return(toString(unique(ann[ann$junction_aa==x,]$Epitope.gene)))}))
exp.b$pep<-unlist(lapply(exp.b$junction_aa, function(x) {return(toString(unique(ann[ann$junction_aa==x,]$Epitope.gene)))}))

exp.a$pep[exp.a$pep==""]<-"n/a"
exp.b$pep[exp.b$pep==""]<-"n/a"


# function to make graph with all annotations I need
prepare_graph<-function(exp.a){
  kma<-kernelMatrix(sk, exp.a$junction_aa, exp.a$junction_aa)
  ones<-which(kma>0.98)
  kma[ones]<-0
  
  connect_a<-(kma>=tha)
  km_adj_a<-kma
  km_adj_a[which(connect_a)]<-1
  km_adj_a[which(!connect_a)]<-0
  
  km_graph_a<-graph.adjacency(km_adj_a, mode="undirected",weighted = TRUE,add.colnames=NULL,diag=FALSE)
  
  V(km_graph_a)$clus<-clusters(km_graph_a)$membership
  V(km_graph_a)$cdr3<-exp.a$junction_aa
  V(km_graph_a)$Vgene<-exp.a$v_call
  V(km_graph_a)$Jgene<-exp.a$j_call
  V(km_graph_a)$ID<-exp.a$ID
  V(km_graph_a)$ann<-exp.a$ann
  V(km_graph_a)$control<-exp.a$control
  V(km_graph_a)$public<-exp.a$public
  if ("pep" %in% names(exp.a)){
    V(km_graph_a)$cd4_cd8<-exp.a$cd4_cd8
    V(km_graph_a)$epitope<-exp.a$pep 
  }
  
  return(km_graph_a)
}


## save function
save_node_info<-function(nodes.a, chain){
  nodes.a.1<-dcast(nodes.a, clus~ID)
  a.1<-nodes.a.1[,2:36]
  nodes.a.1$clus.size<-rowSums(a.1)
  nodes.a.1$num.IDs<-rowSums(a.1!=0)
  date<-Sys.Date()
  # name<- paste(p0, "output/", date, "_", chain, "_cluster_patID.csv", sep="")
  # write.csv(nodes.a.1, name)
  return(nodes.a.1)
}

km_graph_a<-prepare_graph(exp.a)
km_graph_b<-prepare_graph(exp.b)

# km_graph_a.wc<-prepare_graph(exp.a.wc)
# km_graph_b.wc<-prepare_graph(exp.b.wc)

nodes.a<-igraph::as_data_frame(km_graph_a, what="vertices")
nodes.b<-igraph::as_data_frame(km_graph_b, what="vertices")

## investigating clusters that are made of only one ID - alpha

# clusters.a<-names(table(nodes.a$clus)[table(nodes.a$clus) > 2])
# nodes.clus.a<-nodes.a[nodes.a$clus %in% clusters.a,]
# 
# for (c in clusters.a){
#   s<-nodes.clus.a[nodes.clus.a$clus == c,]
#   ids<-unique(s$ID)
#   if (length(ids) == 1){
#     print(paste("cluster", c, "composed of single ID:", ids, ", size =", dim(s)[1], ", unique CDRs:", length(unique(s$cdr3))))
#   }
# }

## investigating clusters that are made of only one ID - beta

# clusters<-names(table(nodes.b$clus)[table(nodes.b$clus) > 2])
# nodes.clus<-nodes.b[nodes.b$clus %in% clusters,]
# 
# for (c in clusters){
#   s<-nodes.clus[nodes.clus$clus == c,]
#   ids<-unique(s$ID)
#   if (length(ids) == 1){
#     print(paste("cluster", c, "composed of single ID:", ids, ", size =", dim(s)[1], ", unique CDRs:", length(unique(s$cdr3))))
#   }
# }
# 
# interesting.ids<-c(90, 153, 175, 267, 268, 324)
# 
# counts<-data.frame(table(nodes.b$ID))
# names(counts)<-c("ID", "cdr3.count")
# 
# counts$interesting.ids<-unlist(lapply(counts$ID,
#                                        function(x) ifelse(x %in% interesting.ids, "yes", "no")))
# 
# ggplot(counts, aes(x = interesting.ids, y = cdr3.count, color = ID, group = interesting.ids)) + 
#   geom_boxplot(width = 0.5, fill = NA) + geom_jitter(width = 0.1, size = 3)
# 
# interesting.clusters<-c(523, 556, 558, 1684, 3191) # only keeping clusters > 4
# 
# mycdr3s<-list()
# 
# for (c in interesting.clusters){
#   s<-nodes.b[nodes.b$clus == c,]
#   cdr3s<-unique(s$cdr3)
#   mycdr3<-list(cdr3s)
#   names(mycdr3)<-c(c)
#   mycdr3s<-c(mycdr3s, mycdr3)
# }

### now check that really the clusters are made of cdr3s from the ID I think

# for (c in interesting.clusters){
#   cdr3s<-mycdr3s[as.character(c)]
#   s<-exp_AB_wide3[exp_AB_wide3$junction_aa %in% unlist(cdr3s),]
#   print(paste("cdr3s in cluster", c, "are from IDs", unique(s$ID)))
# }


# nodes.a.wc<-igraph::as_data_frame(km_graph_a.wc, what="vertices")
# nodes.b.wc<-igraph::as_data_frame(km_graph_b.wc, what="vertices")

nodes.a.1<-save_node_info(nodes.a, "alpha.noctrl")
nodes.b.1<-save_node_info(nodes.b, "beta.noctrl")

# nodes.a.wc.1<-save_node_info(nodes.a.wc, "alpha.wctrl")
# nodes.b.wc.1<-save_node_info(nodes.b.wc, "beta.wctrl")

singlet.doublet.b<-nodes.b.1[nodes.b.1$clus.size<3,]$clus
singlet.doublet.a<-nodes.a.1[nodes.a.1$clus.size<3,]$clus

# singlet.doublet.b.wc<-nodes.b.wc.1[nodes.b.wc.1$clus.size<3,]$clus
# singlet.doublet.a.wc<-nodes.a.wc.1[nodes.a.wc.1$clus.size<3,]$clus

isolates <- function(g){return(which(igraph::degree(g)==0))}
b.singlet.doublet<-function(g){return(which(V(g)$clus %in% singlet.doublet.b))}
a.singlet.doublet<-function(g){return(which(V(g)$clus %in% singlet.doublet.a))}

# b.wc.singlet.doublet<-function(g){return(which(V(g)$clus %in% singlet.doublet.b.wc))}
# a.wc.singlet.doublet<-function(g){return(which(V(g)$clus %in% singlet.doublet.a.wc))}

g.b.s<-delete.vertices(km_graph_b,b.singlet.doublet(km_graph_b))
g.a.s<-delete.vertices(km_graph_a,a.singlet.doublet(km_graph_a))

# g.b.wc.s<-delete.vertices(km_graph_b.wc,b.wc.singlet.doublet(km_graph_b.wc))
# g.a.wc.s<-delete.vertices(km_graph_a.wc,a.wc.singlet.doublet(km_graph_a.wc))

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## plot all kinds of clustering plot I need

plot_all <- function(g.b.s, col_vector, chain){
  e <- get.edgelist(g.b.s,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g.b.s))
  n <- length(unique(V(g.b.s)$Vgene))
  col_vector1<-sample(col_vector,n, replace=F)
  names(col_vector1)<-levels(factor(unique(V(g.b.s)$Vgene)))
  V(g.b.s)$col.V<-unlist(lapply(V(g.b.s), function(x) {as.character(col_vector1[unlist(V(g.b.s)[x]$Vgene)])}))
  V(g.b.s)$shape<-unlist(lapply(V(g.b.s), function(x) {ifelse(V(g.b.s)[x]$control == TRUE, "square", "circle")}))
  shapes<-c("CD4"="blue", "CD8"="red", "CD8, CD4"="green", "na"="darkgrey")
  V(g.b.s)$shape.cd8<-unlist(lapply(V(g.b.s), function(x) as.character(shapes[unlist(V(g.b.s)[x]$cd4_cd8)])))
  
  
  date<-Sys.Date()
#   name<-paste(folder_plots, chain, "_clus_Vgene.svg", sep="")
  
#   svg(name) 
#   par(xpd=F, mar=c(8.5,8.5,8.5,8.5))
#   plot(g.b.s, edge.width=1.5, vertex.label=NA, vertex.size=3.5,
#        layout = l, vertex.color=as.character(vertex_attr(g.b.s, "col.V")),
#        bty="T", vertex.shape=as.character(vertex_attr(g.b.s, "shape"))
#   )
#   legend(1.05,1, legend=names(col_vector1), fill=col_vector1, ncol=2, cex=0.6)
#   dev.off()
  
  ids<-unique(exp_AB_wide3$ID)
  n <- length(ids)
  set.seed(123)
  col_vector2<-sample(col_vector, n, replace=F)
  names(col_vector2)<-levels(factor(ids))
  V(g.b.s)$col.ID<-lapply(V(g.b.s), function(x) {as.character(col_vector2[unlist(as.character(V(g.b.s)[x]$ID))])})
  
  name<-paste(folder_plots, "Fig2f_",  chain, "_clus_ID.svg", sep="")
  svg(name) 
  par(xpd=F, mar=c(8.5,8.5,8.5,8.5))
  plot(g.b.s, edge.width=1.5, vertex.label=NA, vertex.size=3.5,
       layout = l, vertex.color=as.character(vertex_attr(g.b.s, "col.ID")),
       vertex.shape=as.character(vertex_attr(g.b.s, "shape"))
  )
  legend(1.05,1, legend=names(col_vector2), fill=col_vector2, ncol=2, cex=0.6)
  dev.off()
  
#   n <- length(unique(V(g.b.s)$Jgene))
#   col_vector3<-sample(col_vector, n, replace=F)
#   names(col_vector3)<-levels(factor(unique(V(g.b.s)$Jgene)))
#   V(g.b.s)$col.J<-lapply(V(g.b.s), function(x) {as.character(col_vector3[unlist(as.character(V(g.b.s)[x]$Jgene))])})
  
#   name<-paste(p2, date, "_", chain, "clus_Jgene.svg", sep="")
#   svg(name) 
#   par(xpd=F, mar=c(8.5,8.5,8.5,8.5))
#   plot(g.b.s, edge.width=1.5, vertex.label=NA, vertex.size=3.5,
#        layout = l, vertex.color=as.character(vertex_attr(g.b.s, "col.J")),
#        vertex.shape=as.character(vertex_attr(g.b.s, "shape"))
#   )
#   legend(1.05,1, legend=names(col_vector3), fill=col_vector3, ncol=2, cex=0.6)
#   dev.off()
  
#   col_vector4<-c("red", "grey")
#   names(col_vector4)<-c("ann", "not.ann")
#   V(g.b.s)$col.ann<-lapply(V(g.b.s), function(x) {as.character(col_vector4[unlist(as.character(V(g.b.s)[x]$ann))])})
  
#   name<-paste(p2, date, "_", chain, "_clus_ann.svg", sep="")
#   svg(name) 
#   par(xpd=F, mar=c(8.5,8.5,8.5,8.5))
#   plot(g.b.s, edge.width=1.5, vertex.label=NA, vertex.size=3.5,
#        layout = l, vertex.color=as.character(vertex_attr(g.b.s, "col.ann")),
#        vertex.shape=as.character(vertex_attr(g.b.s, "shape"))
#   )
#   legend(1.05,1, legend=names(col_vector4), fill=col_vector4, ncol=2, cex=0.6)
#   dev.off()
  
#   my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 10000)
  
#   sf <- max(abs(V(g.b.s)$cd8))
#   node.colors <- (V(g.b.s)$cd8+sf) / (2*sf) * 10000
  
#   n <- length(unique(V(g.b.s)$public))
#   col_vector6<-gg_color_hue(n)
#   names(col_vector6)<-levels(factor(unique(V(g.b.s)$public)))
#   col_vector6[1]<-"white"
#   V(g.b.s)$col.pub<-lapply(V(g.b.s), function(x) {as.character(col_vector6[unlist(as.character(V(g.b.s)[x]$public))])})
  
#   name<-paste(folder_plots, chain, "_clus_publicity.svg", sep="")
#   svg(name) 
#   par(xpd=F, mar=c(8.5,8.5,8.5,8.5))
#   plot(g.b.s, edge.width=1.5, vertex.label=NA, vertex.size=3.5,
#        layout = l, vertex.color=as.character(vertex_attr(g.b.s, "col.pub")),
#   )
#   legend(1.05,1, legend=names(col_vector6), fill=col_vector6, ncol=2, cex=0.6)
#   dev.off()
  
  eps<-unique(ann$Epitope.gene)
  
  if ("epitope" %in% names(vertex_attr(g.b.s))){
    n <- length(unique(eps))
    col_vector5<-gg_color_hue(n)
    names(col_vector5)<-levels(factor(unique(eps)))
    col_vector5["n/a"]<-"grey"
    V(g.b.s)$col.ep<-lapply(V(g.b.s)$epitope, function(x) {as.character(col_vector5[unlist(as.character(x))])})
    num.eps<-unlist(lapply(V(g.b.s)$epitope, function(x)length(unlist(strsplit(as.character(x), split=", ")))))
    assert("no more than 2 epitopes assigned", max(unlist(num.eps)) == 2)
    values<-lapply(num.eps, function(x) return(c(0,0)))
    values1<-lapply(1:length(num.eps), function(x) {unlist(ifelse(num.eps[x]==1, list(c(1,0)), list(c(1, 1))))})
    cols<-lapply(V(g.b.s)$epitope, function(x) {unlist(ifelse(length(unlist(strsplit(x, split=", "))) > 1, 
          list(c(as.character(col_vector5[unlist(strsplit(x, split=", "))[1]]), as.character(col_vector5[unlist(strsplit(x, split=", "))[2]]))),
          list(c(as.character(col_vector5[unlist(strsplit(x, split=", "))[1]]), "black"))))
    })
    V(g.b.s)$pie.color<-cols
    V(g.b.s)$pie.values<-values1
      
    name<-paste(folder_plots, "FigS13_", chain, "_clus_epitope.svg", sep="")
    svg(name) 
    par(xpd=F, mar=c(8.5,8.5,8.5,8.5))
    plot(g.b.s, edge.width=1.5, vertex.label=NA, vertex.size=3.5,
         layout = l, 
         vertex.shape="pie", vertex.pie= values1
    )
    legend(1.05,1, legend=names(col_vector5), fill=col_vector5, ncol=2, cex=0.6)
    dev.off()
  }
  
  return(g.b.s)
}

g.a.s<-plot_all(g.a.s, col_vector, "alpha.noctrl_3orMore")
g.b.s<-plot_all(g.b.s, col_vector, "beta.noctrl_3orMore")

####

# print sequence logos

# g.a.17<-igraph::induced_subgraph(g.a.s, V(g.a.s)[V(g.a.s)$clus==17])
# 
# plot(g.a.17, edge.width=2, vertex.label=NA, vertex.size=7,
#      layout = layout.kamada.kawai,
#      vertex.shape="pie", vertex.pie=vertex_attr(g.a.17, "pie.values"))
# 
# align_a17 <- msa(unlist(unique(V(g.a.17)$cdr3)),type="protein")
# conMat_as_u<-consensusMatrix(align_a17)
# ggseqlogo(conMat_as_u, method="probability")
# 
# align_a17.1 <- msa(unlist(V(g.a.17)$cdr3),type="protein")
# conMat_as<-consensusMatrix(align_a17.1)
# ggseqlogo(conMat_as, method="probability")
# 
# a17<-igraph::as_data_frame(g.a.17, what="vertices")
# View(a17)
# 
# dcast(a17, ID~c(ann,cd4_cd8, epitope))

#####
# 
# ann1<-dcast(ann, junction_aa~Epitope.gene)
# ann1[2:dim(ann1)[2]][ann1[2:dim(ann1)[2]] > 1] <-1
# ann1$num.eps<-rowSums(ann1[2:dim(ann1)[2]])
# 
# cdr3.mul<-ann1[ann1$num.eps>1,]$junction_aa
# 
# repeats<-ann[ann$junction_aa %in% cdr3.mul,]
# write.csv(repeats, paste(folder_data, "exp_AB_ann_multiple_epitopes.csv", sep=""))
# 
# # fiveormore.b<-nodes.b.1[nodes.b.1$clus.size<5,]$clus
# # fivers.b<-function(g){return(which(V(g)$clus %in% fiveormore.b))}
# # g.b.s1<-delete.vertices(km_graph_b,fivers.b(km_graph_b))
# # 
# # fiveormore.a<-nodes.a.1[nodes.a.1$clus.size<5,]$clus
# # fivers.a<-function(g){return(which(V(g)$clus %in% fiveormore.a))}
# # g.a.s1<-delete.vertices(km_graph_a,fivers.a(km_graph_a))
# # 
# # plot_all(g.a.s1, col_vector, "alpha_5orMore")
# # plot_all(g.b.s1, col_vector, "beta_5orMore")

####
# look at clustering in non-expanded (makes equivalent cluster networks)

# load(paste(p1, "control_alpha_1.RData", sep=""))
# load(paste(p1, "control_beta_1.RData", sep=""))
# 
# prepare_graph_red<-function(exp.a){
#   kma<-kernelMatrix(sk, exp.a$junction_aa, exp.a$junction_aa)
#   ones<-which(kma>0.98)
#   kma[ones]<-0
#   
#   connect_a<-(kma>=tha)
#   km_adj_a<-kma
#   km_adj_a[which(connect_a)]<-1
#   km_adj_a[which(!connect_a)]<-0
#   
#   km_graph_a<-graph.adjacency(km_adj_a, mode="undirected",weighted = TRUE,add.colnames=NULL,diag=FALSE)
#   
#   V(km_graph_a)$clus<-clusters(km_graph_a)$membership
#   V(km_graph_a)$cdr3<-exp.a$junction_aa
#   V(km_graph_a)$ann<-exp.a$ann
#   
#   return(km_graph_a)
# }
# 
# save_node_info1<-function(nodes.a, chain){
#   nodes.a.1<-dcast(nodes.a, clus~ann)
#   a.1<-data.frame(nodes.a.1[,2:dim(nodes.a.1)[2]])
#   # print(a.1)
#   nodes.a.1$clus.size<-rowSums(a.1)
#   nodes.a.1$num.IDs<-rowSums(a.1!=0)
#   date<-Sys.Date()
#   # name<- paste(p0, "output/", date, "_", chain, "_cluster_ann.csv", sep="")
#   # write.csv(nodes.a.1, name)
#   return(nodes.a.1)
# }
# 
# plot_red <- function(g.b.s, col_vector, chain){
#   e <- get.edgelist(g.b.s,names=FALSE)
#   l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g.b.s))
# 
#   date<-Sys.Date()
#   
#   col_vector4<-c("red", "grey")
#   names(col_vector4)<-c("ann", "not.ann")
#   V(g.b.s)$col.ann<-lapply(V(g.b.s), function(x) {as.character(col_vector4[unlist(as.character(V(g.b.s)[x]$ann))])})
#   
#   name<-paste(folder_plots, chain, "clus_ann.svg", sep="")
#   # svg(name) 
#   par(xpd=F, mar=c(8.5,8.5,8.5,8.5))
#   plot(g.b.s, edge.width=1.5, vertex.label=NA, vertex.size=3.5,
#        layout = l, vertex.color=as.character(vertex_attr(g.b.s, "col.ann")),
#   )
#   legend(1.05,1, legend=names(col_vector4), fill=col_vector4, ncol=2, cex=0.6)
#   # dev.off()
# }
# 
# 
# for (i in 1:length(control_a)){
#   print(i)
#   ne.a<-control_a[i]
#   ne.a<-data.frame(ne.a)
#   names(ne.a)<-c("junction_aa")
#   ne.a$ann<-unlist(lapply(ne.a$junction_aa, function(x) ifelse(x %in% all.ann.a, "ann", "not.ann")))
#   
#   km_graph_a.ne<-prepare_graph_red(ne.a)
#   nodes.a.ne<-igraph::as_data_frame(km_graph_a.ne, what="vertices")
#   nodes.a.ne.1<-save_node_info1(nodes.a.ne, paste("alpha.nonexp", i, sep=""))
#   singlet.doublet.a.ne<-nodes.a.ne.1[nodes.a.ne.1$clus.size<3,]$clus
#   a.ne.singlet.doublet<-function(g){return(which(V(g)$clus %in% singlet.doublet.a.ne))}
#   
#   g.a.ne.s<-delete.vertices(km_graph_a.ne,a.ne.singlet.doublet(km_graph_a.ne))
#   print("plotting...")
#   if (length(V(g.a.ne.s))>0){
#     print("plotting...")
#     plot_red(g.a.ne.s, col_vector, paste("alpha.nonexp", i, sep="")) 
#   }
# }
# 
# for (i in 1:length(control_b)){
#   print(i)
#   ne.b<-control_b[i]
#   ne.b<-data.frame(ne.b)
#   names(ne.b)<-c("junction_aa")
#   ne.b$ann<-unlist(lapply(ne.b$junction_aa, function(x) ifelse(x %in% all.ann.b, "ann", "not.ann")))
#   
#   print("graphing...")
#   km_graph_b.ne<-prepare_graph_red(ne.b)
#   print("save nodes...")
#   nodes.b.ne<-igraph::as_data_frame(km_graph_b.ne, what="vertices")
#   print("save csv...")
#   nodes.b.ne.1<-save_node_info1(nodes.b.ne, paste("beta.nonexp", i, sep=""))
#   singlet.doublet.b.ne<-nodes.b.ne.1[nodes.b.ne.1$clus.size<3,]$clus
#   b.ne.singlet.doublet<-function(g){return(which(V(g)$clus %in% singlet.doublet.b.ne))}
#   
#   print("remove singlets and doublets")
#   g.b.ne.s<-delete.vertices(km_graph_b.ne,b.ne.singlet.doublet(km_graph_b.ne))
#   # print(g.b.ne.s)
#   
#   if (length(V(g.b.ne.s))>0){
#     print("plotting...")
#     plot_red(g.b.ne.s, col_vector, paste("beta.nonexp", i, sep="")) 
#   }
# }
