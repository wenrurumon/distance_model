setwd('/home/zhu/rushdata/methylation_net')
load('genesel_rawi.rda')
library(fda)
library(flare)

ftransform <- function(X,pos,nbasis=31){
	pos <- (pos-min(pos))/(max(pos)-min(pos))
	X <- X[,order(pos),drop=F]
	pos <- pos[order(pos)]
	if(is.na(pos)){
		X <- cbind(X,X);
		pos <- 0:1
	}
	fbasis<-create.fourier.basis(c(0,1))
	fphi <- eval.basis(pos,fbasis)
	fcoef <- ginv(t(fphi)%*%fphi)%*%t(fphi)%*%t(X)
	rlt <- t(fcoef-rowMeans(fcoef))/sqrt(nrow(X))
	return(rlt)
}

model <- function(i){
	X <- genesel_rawi[[i]][[1]]
	pos <- X <- genesel_rawi[[i]][[2]]
	ftransform(X,pos,31)
}
