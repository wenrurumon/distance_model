
rm(list=ls())
setwd('/home/zhu/rushdata/methylation_net')
load('snpbycol.rda')

gene.dist <- dist(t(gene.out))
gene.hc <- hclust(gene.dist)
save(gene.dist,gene.hc,gene.out,file='temp.rda')
