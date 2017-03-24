rm(list=ls())
setwd('/home/zhu/rushdata/methylation_net')
load('genesel_rawi.rda')
load('/home/zhu/rushdata/qtl20170226_paths_methy_snp/pid.rda')
library(fda)
library(flare)

gname <- names(gene.sel)
gene.sel <- lapply(gene.sel,function(x){
	snp <- x[[1]]
	snp <- snp[match(pid,rownames(snp)),,drop=F]
	return(list(snp=snp,pos=x[[2]]))
})

ftransform <- function(X,pos,nbasis=31){
	pos <- (pos-min(pos))/(max(pos)-min(pos))
	X <- X[,order(pos),drop=F]
	pos <- pos[order(pos)]
	if(length(pos)==1){
		X <- cbind(X,X);
		pos <- 0:1
	}
	fbasis<-create.fourier.basis(c(0,1),nbasis)
	fphi <- eval.basis(pos,fbasis)
	fcoef <- ginv(t(fphi)%*%fphi)%*%t(fphi)%*%t(X)
	rlt <- t(fcoef-rowMeans(fcoef))/sqrt(nrow(X))
	return(rlt)
}

model <- function(i,nbasis){
	X <- gene.sel[[i]][[1]]
	pos <- gene.sel[[i]][[2]]
	ftransform(X,pos,nbasis)
}

qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}

gene_51basis <- lapply(1:length(gene.sel),model,nbasis=51)
gene.pool <- do.call(rbind,gene_51basis)
gene.qpca1 <- qpca(gene.pool)
gene.qpca2 <- qpca(gene.pool,rank=which(gene.qpca1$prop>0.9)[1])
gene.qpca2 <- gene.qpca2$X[,1:which(gene.qpca2$prop>0.9)[1],drop=F]

map <- rep(1:length(gene_51basis),each=nrow(gene_51basis[[1]]))
gene.out <- lapply(unique(map),function(i){
	gene.qpca2[map==i,]
})
names(gene.out) <- gname
gene.out <- do.call(cbind,lapply(gene.out,as.vector))

