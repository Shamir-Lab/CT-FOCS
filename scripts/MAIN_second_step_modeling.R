##################################################################################################################################
################################# Second step modeling - inferring cell type-specific EP links ###################################
##################################################################################################################################


########################################## Note: Only works via R unix platform ##################################################

## This is a line by line running tutorial. In order to to run a subset of lines please select them using shift and run using ctrl+R

## Please follow the instruction given in http://acgt.cs.tau.ac.il/ct-focs/tutorial.html

# load libraries - make sure you pre-installed these packages
cran_libs = c('nlme','MASS','parallel','glmmTMB','RColorBrewer','ggplot2')
bioconductor_libs = c('GenomicRanges','biomaRt','topGO')
for (p in (c(bioconductor_libs,cran_libs))){library(p,character.only=T)}


# Variables
k <- 10 #k closest enhancers to each gene
win.size <- 5*10**5 #maximum window size upstream/downstream to select genes with at least 10 enhancers
# Number of parallel cores to run on - use detectCores(all.tests = FALSE, logical = TRUE) to identify the number of cores
mc.cores <- 40
# FDR correction method
method='BH'
# FDR threshold
ranef.pv.thr=0.1


# load objects and read files
data_directory = 'data/'
script_directory = 'scripts/'
tmp_directory = 'tmp/'
# data types - fantom5, encode, roadmap
data_type <- 'encode'

#load code
source(paste(script_directory,'FUNCTIONS_reg_functions.R',sep=''))


# Enhancer genomic positions
enh.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.enh.pos.rds',sep=''))
# Promoters positions
prom.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.prom.pos.rds',sep=''))
# Enhancer count matrix:
Me = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.enh.count.rds',sep='')))
# Promoter count matrix:
Mg = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.prom.count.rds',sep='')))
offs = NULL

if(data_type=='fantom5' | data_type=='roadmap'){
	# Enhancer RLE matrix:
	Me = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.enh.rle.rds',sep='')))
	# Promoter RLE matrix:
	Mg = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.prom.rle.rds',sep='')))
}

# Sample annotations
sample_annot <- readRDS(paste(data_directory,data_type,'/',data_type,'.sample.annot.rds',sep=''))
# Promoter to closest enhnacers in fixed window:
g_to_e = readRDS(paste(data_directory,data_type,'/',data_type,'.cand.enh.rds',sep=''))

#Load the list of robust promoter models from the first step modeling
gene_val_groups <- readRDS(file=paste(tmp_directory,data_type,'.gene_val_groups.rds',sep=''))
#we select only promoter models that passed either both binary and activity level tests or only the activity level test
length(gene_list <- union(gene_val_groups[["Both"]],gene_val_groups[["Level only"]]))

if(data_type=='encode'){
	offs = sample_annot$eff.lib.size
	names(offs)=rownames(sample_annot)
}

# Samples (columns in matrices) library size:
col2type <- sample_annot[,'col2type',drop=F]
col2type[,1] <- as.character(col2type[,1])


tmpMg = Mg
tmpMe = Me

Mg <- as.matrix(tmpMg)
Me <- as.matrix(tmpMe)
Mg <- Mg[,rownames(col2type)]
Me <- Me[,rownames(col2type)]


if(data_type=='fantom5'){
	Mg <- t(t(Mg)+2e6/sample_annot[colnames(tmpMg),]$lib_size)
	Me <- t(t(Me)+2e6/sample_annot[colnames(tmpMg),]$lib_size)
	Mg <- as.matrix(log2(Mg))
	Me <- as.matrix(log2(Me))
}

Mg_fpkm <- Mg

#mark for each gene its intronic TREs
Mg_normalized = Mg
Me_normalized = Me

c2t <- col2type[,1]


## FANTOM5 dataset is log-normal-like distributed, therefore, we use the regular lmm with errors/random_effects assumed normal distributed

if(data_type=='fantom5' | data_type=='roadmap'){
	
	lmm_results <- mclapply(1:length(gene_list),function(i) get_lmm_results(i,gene_list[i],Mg[gene_list[i],],
						Mg,Me,g_to_e,col2type[,1],ranef.pv.thr=ranef.pv.thr,k=k,threshold=0,loop=T,only.pos=T,chau.pv.thr=0.5,MADCutOff=2.5,offs=offs),mc.cores=20)
	
	names(lmm_results) <- gene_list
	gc()
	
	
	saveRDS(lmm_results,file=paste0(tmp_directory,data_type,'.mixed_effects_results.rds'))

}


## ENCODE dataset is Zero-Inflated negative binomial distributed, therefore, we use the glmmTMB with ZINB as the link function

if(data_type=='encode'){
	chunks <- split(gene_list, ceiling(seq_along(gene_list)/500))
	length(n <- 1:length(chunks))
	
	for( j in n){
		print(paste0("Now working on chunk: ",j))
		if(file.exists(paste0(tmp_directory,data_type,'.glmm_res.chunk.',j,'.rds'))) next;
		chunk <- chunks[[j]]
		glmm_res <- mclapply(1:length(chunk),function(i) get_glmm_nb_results(i,chunk[i],Mg_fpkm[chunk[i],],
	                                                                        Mg.count,Me.count,g_to_e,c2t,k=k,threshold=0,offs=offs,MADCutOff=2.5),mc.cores=40)
		names(glmm_res) <- chunk
		saveRDS(glmm_res,file=paste0(tmp_directory,data_type,'.glmm_res.chunk.',j,'.rds'))
		rm(glmm_res)
		gc()
	}
	
	glmm_results <- list()
	for( j in n){
		glmm_results <- c(glmm_results,readRDS(file=paste0(tmp_directory,data_type,'.glmm_res.chunk.',j,'.rds')))
	}
	
	saveRDS(glmm_results,file=paste0(tmp_directory,data_type,'.mixed_effects_results.rds'))
	
	gc()

}
## Analyze ENCODE/FANTOM5 glmm results
## Here glmm_results stands also for the lmm_results

glmm_results = readRDS(file=paste0(tmp_directory,data_type,'.mixed_effects_results.rds'))

glmm_results <- mclapply(glmm_results,function(x) if(x$success){ return(x) }else{ return(NULL)}, mc.cores=10)
glmm_results[sapply(glmm_results,is.null)] = NULL

int.p.val.glmm <- do.call(c,lapply(glmm_results,function(x) if(x$success) return(x$p.v.e[1]) else return(NA)))

summary(p.adjust(int.p.val.glmm,method=method))
table(p.adjust(int.p.val.glmm,method=method)<=ranef.pv.thr)

slope.p.val.glmm <- do.call(c,lapply(glmm_results,function(x) if(x$success) return(x$p.v.e[-1]) else return(NA)))

summary(p.adjust(slope.p.val.glmm,method=method))
table(p.adjust(slope.p.val.glmm,method=method)<=ranef.pv.thr)

# number of tests conducted ~|P|*(k+1) - |P| number of promoters, k=10 closest enhancers + intercept
n = sum(length(slope.p.val.glmm) + length(int.p.val.glmm), 0 - (table(is.na(int.p.val.glmm))[2]+table(is.na(slope.p.val.glmm))[2]),na.rm=T)

# find active P cell types
p.active <- lapply(seq_len(length(glmm_results)),function(i){x<-glmm_results[[i]]; is.sig=TRUE; p.v.e<-p.adjust(x$p.v.e,n=n,method=method); if(is.na(p.v.e[1]) | p.v.e[1]>ranef.pv.thr){is.sig=FALSE};
if(length(x$mad_res$intercept)==0){is.sig=FALSE}; return(list(is.sig=is.sig,cells=x$mad_res$intercept))})
names(p.active) <- names(glmm_results)

# find active E cell types
e.active <- lapply(seq_len(length(glmm_results)),function(i){x<-glmm_results[[i]]; res <- x$mad_res[-1]; p.v.e<-p.adjust(x$p.v.e,n=n,method=method)[-1]; 
		ind=(p.v.e<=ranef.pv.thr); ind[is.na(ind)]=FALSE; return(res[ind])})
names(e.active) <- names(glmm_results)


# create cell to active P list
p2c.active <- sapply(seq_len(length(p.active)),function(i){x<-p.active[[i]]; if(!x$is.sig) return(NULL); if(length(x$cells)==0) return(NULL); x$cells},simplify=F)
names(p2c.active) <- names(p.active)
p2c.active[sapply(p2c.active, is.null)] <- NULL
c2p.active <- inverseList(p2c.active)
e2c.active = as.list(unlist(e.active))
c2e.active <- inverseList(e2c.active)


num.p.active <- sapply(c2p.active,length)
#show how many promoters we have per cell type
sort(num.p.active)
#show how many linked enhancers we have per cell type
sort(sapply(c2e.active,length))


# find for specific cell type its active P and E
shift.size <- (end(prom.bs)-start(prom.bs))/2
prom.bs.mid <- GenomicRanges::shift(prom.bs,shift=shift.size)
strand(prom.bs.mid) <- '+'

shift.size <- (end(enh.bs)-start(enh.bs))/2
enh.bs.mid <- GenomicRanges::shift(enh.bs,shift=shift.size)
#enh.bs.mid <-  GenomicRanges::promoters(enh.bs.shift, upstream=0,downstream=1)
strand(enh.bs.mid) <- '+'


cells <- names(c2p.active)
df.list.focs <- data.frame()


res <- mclapply(seq_len(length(cells)),function(i) createDfList(i,cells[i],c2p.active[[cells[i]]],e.active), mc.cores=20)

df.list.focs = do.call(rbind,lapply(res,function(x) x[[1]]))
ep_links = lapply(res,function(x) x[[2]])
names(ep_links) = cells

ep_links <- sapply(names(ep_links),function(x) {y=ep_links[[x]]; y[sapply(y,is.null)] <- NULL; return(y)})

# how many EP links we have per cell type
sort(sapply(ep_links,function(x) length(unlist(x))),decreasing=F)

saveRDS(ep_links,file=paste0(data_directory,data_type,'.EP_links.cell.rds'))

### find EP links for a specific cell type

ep_links = readRDS(file=paste0(data_directory,data_type,'.EP_links.rds'))

#what cell types do we have?
(cells <- names(ep_links))

s.cell <- 'GM12878'

if(data_type=='fantom5')
	s.cell = tolower(s.cell)

e.active.cell <- ep_links[[s.cell]]
p.active.cell <- names(e.active.cell)
length(p.active.cell)

tmp <- e.active.cell
tmp[sapply(tmp,is.null)] <- NULL

length(unlist(tmp)) # num of EP links in s.cell
length(unique(unlist(tmp))) # num of linked enhancers in s.cell
length(tmp) # num of linked promoters in s.cell


