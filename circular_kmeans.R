
rm(list=ls())
load("C:/Users/zhu2/Application Data/SSH/temp/temp (4).rda")

dist2 <- as.matrix(gene.dist)^2
kmeans_initn <- function(dist2){
  kmi <- kmeans(dist2,2,iter.max=10)
  grp <- list(dist2[kmi$cluster==1,kmi$cluster==1],
              dist2[kmi$cluster!=1,kmi$cluster!=1])
  return(grp)
}

rlti <- kmeans_initn(dist2)
while(max(sapply(rlti,ncol))>80){
  i <- which(sapply(rlti,ncol)==max(sapply(rlti,ncol)))[1]
  disti <- rlti[[i]]
  rlti <- rlti[-i]
  rlti <- c(rlti,kmeans_initn(disti))
  print(sapply(rlti,ncol))
}
