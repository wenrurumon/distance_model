rm(list=ls())
data <- lapply(1:22,function(args){
	print(args)
	datafolder <- '/home/zhu/rushdata/rawdata/methylation_rawdata/methylation_imputed_all_chr/'
	mapfolder <- '/home/zhu/rushdata/rawdata/methylation_rawdata/illqc_probes_loc/'
	datafile <- paste0(datafolder,'dnaMeth_matrix_chr',args,'_748qc_imputed.txt')
	posfile <- paste0(mapfolder,'illqc_probes_chr',args,'_loc.txt')
	mapfile <- '/home/zhu/rushdata/rawdata/methylation_rawdata/GPL13534.map'
	idfile <- '/home/zhu/rushdata/rawdata/methylation_rawdata/id_delt.txt'
	raw <- t(read.table(datafile))
	pos <- read.table(posfile,header=T)
	map <- read.table(mapfile)
	id <- read.table(idfile,header=T)
	rownames(raw) <- id$projid
	list(raw=raw,pos=pos)
	})

setwd('/home/zhu/rushdata/code')
sel <- readLines('methylation_candidate.csv')
pos <- lapply(data,function(x) x[[2]])
map <- read.table('/home/zhu/rushdata/rawdata/methylation_rawdata/GPL13534.map')
map.sel <- map[map$V5%in%sel,]

data.m <- do.call(cbind,lapply(data,function(x) x[[1]]))
data.p <- do.call(rbind,lapply(data,function(x) x[[2]]))

gene <- paste(unique(map.sel$V5))

gene.sel <- lapply(gene,function(genei){
	print(genei)
	pos.sel <- which(paste(data.p$TargetID)%in%paste(map.sel$V4[map.sel$V5==genei]))
	rawi <- data.m[,pos.sel,F]
	posi <- data.p[pos.sel,,F]
	return(list(rawi,posi))
})
setwd('/home/zhu/rushdata/methylation_net')
save(gene.sel,file='genesel_rawi.rda')

############################

rm(list=ls())
setwd('/home/zhu/rushdata/processed/methylation')
f <- dir()
pmethy <- lapply(f,function(x){
	load(x)
	pmethy
	})
pmethy <- do.call(c,pmethy)
setwd('/home/zhu/rushdata/code')
sel <- readLines('methylation_candidate.csv')
pmethy <- pmethy[names(pmethy)%in%sel]
pmethy <- lapply(pmethy,function(x){x[1:2]})

load('/home/zhu/rushdata/methylation_net/genesel_rawi.rda')
