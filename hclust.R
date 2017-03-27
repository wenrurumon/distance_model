rm(list=ls())
setwd('/home/zhu/rushdata/methylation_net')
load("temp.rda")

cluster <- cutree(gene.hc,500)
cluster <- lapply(unique(cluster),function(i){
  names(which(cluster==i))
 })

data.cluster <- lapply(cluster,function(i){
  gene.out[,colnames(gene.out)%in%i,drop=F]
})
