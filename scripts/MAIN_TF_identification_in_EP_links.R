##################################################################################################################################
################################# Identification of known TFs in cell type-specific EP links #####################################
##################################################################################################################################

########################################## Note: Only works via R unix platform ##################################################

## This is a line by line running tutorial. In order to to run a subset of lines please select them using shift and run using ctrl+R

## Please follow the instruction given in http://acgt.cs.tau.ac.il/ct-focs/tutorial.html

## Code for figures - available at MAIN_paper_figures.R script

########################################### Importatnt - MEME suit installation ##################################################
## you should have an installed MEME suit programs under your Linux OS: http://meme-suite.org/doc/download.html
## After installing the MEME suit you should download the MEME's motif databases: http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz
## The motif databases folder (motif_databases) should be located within meme/db/ folder
## Download the file 'HOCOMOCOv11_core_annotation_HUMAN_mono.tsv' from CT-FOCS download page under 'Additional data' table. Insert this file to meme/db/motif_databases/HUMAN/ folder.
##################################################################################################################################

# load libraries - make sure you pre-installed these packages
cran_libs = c('S4Vectors','MASS','parallel','seqinr','MatchIt','RColorBrewer','ggplot2','gplots')
bioconductor_libs = c('GenomicRanges','BSgenome.Hsapiens.UCSC.hg19')
for (p in (c(bioconductor_libs,cran_libs))){library(p,character.only=T)}


# Variables
k <- 10 #k closest enhancers to each gene
win.size <- 0.5e3 #size of upstream/downstream flank size (window = 1000 bp) for the loop's anchors and EP regions
# Number of parallel cores to run on - use detectCores(all.tests = FALSE, logical = TRUE) to identify the number of cores
mc.cores <- 40

# load objects and read files
data_directory = 'data/'
script_directory = 'scripts/'
tmp_directory = 'tmp/'
# data types - fantom5, encode
data_type <- 'encode'

## Motif finding variables
# MEME path
meme_path <- '../meme/'
pwm_path <- paste0(meme_path,'db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme')
ratio <- 4 # number of control matched sequences to each target sequence
pv.th = 0.05 # HyperGeometric p value threshold of significant enriched transcription factor (TF)

if(data_type=='fantom5')
	ratio <- 2

EP_num_th <- 50 #Minimum number of cell type-specific EP links for identifying cell type-specific TFs 
genome <- BSgenome.Hsapiens.UCSC.hg19 # for fetching FASTA sequences based on genomic positions in GRanges object

## Motif finding source functions
source(paste0(script_directory,'FUNCTIONS_motif_finding.R'))

## Add your meme bin path to your PATH
Sys.setenv(PATH = paste(Sys.getenv("PATH"),paste0(meme_path,'bin'),sep = ":"))

## Motif finding functions
source(paste(script_directory,'FUNCTIONS_motif_finding.R',sep=''))

k <- 10 #k closest enhancers to each gene
win.size <- 5*10**5 #maximum window size upstream/downstream to select genes with at least 10 enhancers
mc.cores <- 40
ranef.pv.thr=0.1

## Load 8.4M digital genomic footprints (hg19 genome build) from Neph et al., Nature (2012)
dgf.bs <- readRDS(file=paste0(data_directory,"dgf.bs.rds"))
# Enhancer genomic positions
enh.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.enh.pos.rds',sep=''))
# Promoters positions
prom.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.prom.pos.rds',sep=''))
# Sample annotations
sample.annot <- readRDS(paste(data_directory,data_type,'/',data_type,'.sample.annot.rds',sep=''))
# Promoter to closest enhnacers in fixed window:
g_to_e = readRDS(paste(data_directory,data_type,'/',data_type,'.cand.enh.rds',sep=''))
# Predicted cell type-specific EP links from CT-FOCS
ep_links = readRDS(file=paste0(data_directory,data_type,'/',data_type,'.EP_links.cell.rds'))
#Load the list of robust promoter models from the first step modeling
gene_val_groups <- readRDS(file=paste(tmp_directory,data_type,'.gene_val_groups.rds',sep=''))
#we select only promoter models that passed either both binary and activity level tests or only the activity level test
length(gene_list <- union(gene_val_groups[["Both"]],gene_val_groups[["Level only"]]))

## Some of the TF names in the PWM file are wrong, therefore, we correct the names using the HOCOMOCOv111 TF annotation file
df = read.delim(paste0(data_directory,meme_path,'db/motif_databases/HUMAN/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv'), header=TRUE)
tf_hoc_to_tf = as.character(df[,1])
tf_hoc_to_tf = c(tf_hoc_to_tf,"EVX2_HUMAN.H11MO.0.A") ## add some missing TFs in the file
tf_hoc_to_tf = sapply(tf_hoc_to_tf,function(x) unlist(strsplit(x,"_"))[1])
tf.names <- as.character(df$Transcription.factor)
tf.names <- c(tf.names,'EVX2')
tf.widths <- df$Model.length
tf.widths <- c(tf.widths,13)
names(tf.widths) <- tf.names

df.annot.hoc = data.frame(TF=tf.names,Length=tf.widths)
rownames(df.annot.hoc) = tf_hoc_to_tf

## Create the set of all candidate EP links from the robust promoter models achieved after the first step of modeling
res <- lapply(1:length(gene_list), function(i) g_to_e[[gene_list[i]]][1:k])
names(res) <- gene_list
ep_links_all <- res

prom.bs.all <- prom.bs[names(ep_links_all)]
enh.bs.all <- enh.bs[unique(unlist(ep_links_all))]

# create E/P sequences
seqs.bg.p <- as.list(getSeq(genome, prom.bs.all, as.character=TRUE))
seqs.bg.e <- as.list(getSeq(genome, enh.bs.all, as.character=TRUE))

# create X - the covariate matrix
# X contains one row per enhancer sequence
# Each row contain the counts of A,C,G,T and all possible pairs (e.g., AT,AC,CG,GC...)
X.e <- createMatchItCovMatrix(seqs.bg.e)

# X contains one row per promoter sequence
X.p <- createMatchItCovMatrix(seqs.bg.p)

## Find DGFs within enahncers and promoters in EP links
hit.p <- findOverlaps(dgf.bs,prom.bs.all,ignore.strand=TRUE,type="within")
hit.e <- findOverlaps(dgf.bs,enh.bs.all,ignore.strand=TRUE,type="within")
length(dgf.id <- union(queryHits(hit.p),queryHits(hit.e)))
dgf.bs.reduced.all <- reduce(dgf.bs[sort(dgf.id)])
dgf.bs.reduced.all.p <- reduce(dgf.bs[sort(unique(queryHits(hit.p)))])
dgf.bs.reduced.all.e <- reduce(dgf.bs[sort(unique(queryHits(hit.e)))])

						
cells=names(which(sort(sapply(ep_links,function(x) length(unlist(x))),decreasing=T)>=EP_num_th))

list.res = rep(list(NULL),length(cells))
names(list.res) = cells

list.tf.count = rep(list(NULL),length(cells))
names(list.tf.count) = cells

list.tf.eF = rep(list(NULL),length(cells))
names(list.tf.eF) = cells

list.tf.pVal = rep(list(NULL),length(cells))
names(list.tf.pVal) = cells

######################################################## FIMO analysis #######################################################
last.pos = 1
for(i in 1:length(cells)){
	set.seed(1)
	s.cell <- cells[i]
	print(paste0(i, " Now working on cell: ", s.cell))
	
	e.active.cell <- ep_links[[s.cell]]
	p.active.cell <- names(e.active.cell)
	length(p.active.cell)
	
	
	tmp <- e.active.cell
	tmp[sapply(tmp,is.null)] <- NULL
	length(unlist(tmp)) # num of EP
	length(unique(unlist(tmp))) # num of E
	length(tmp) # num of P
	
	# get the reduced promoter genomic positions according to ep_links
	prom.bs.ep <- prom.bs[names(tmp)]
	enh.bs.ep <- enh.bs[unique(unlist(tmp))]
	
	# find overlaps with dgf	
	hit.p <- findOverlaps(dgf.bs,prom.bs.ep,ignore.strand=TRUE,type="within")
	hit.e <- findOverlaps(dgf.bs,enh.bs.ep,ignore.strand=TRUE,type="within")

	# how many DGFs fall within the enhancers and promoters
	length(dgf.id <- union(queryHits(hit.p),queryHits(hit.e)))
	dgf.bs.reduced.p <- dgf.bs[sort(queryHits(hit.p))]
	dgf.bs.reduced.e <- dgf.bs[sort(queryHits(hit.e))]
	
	length(dgf.bs.reduced.p)
	length(dgf.bs.reduced.e)
	
	
	# run FIMO on target sequences
	fimo.gr.e <- runFimo(dgf.bs.reduced.e,genome=genome, type=data_type)
	fimo.gr.p <- runFimo(dgf.bs.reduced.p,genome=genome, type=data_type)
	
	# Tr is a vector marking the target enhancers (labeled as 1) and the control enhancers (labeled as 0)
	Tr.e <- createMatchItTargetVector(X.e, seqs.bg.e, enh.bs.ep)
	
	# and for promoters
	Tr.p <- createMatchItTargetVector(X.p, seqs.bg.p, prom.bs.ep)
	
	
	#enhancers
	set.id.e.crt <- createMatchItBgSet(X.e, Tr.e, ratio)
	
	#promoters
	set.id.p.crt <- createMatchItBgSet(X.p, Tr.p, ratio)
	
	
	# get the reduced promoter genomic positions according to ep_links
	prom.bs.ep.crt <- prom.bs[names(seqs.bg.p[set.id.p.crt])]
	enh.bs.ep.crt <- enh.bs[names(seqs.bg.e[set.id.e.crt])]
	prom.bs.ep.crt$is_tg = rep(0, length(prom.bs.ep.crt))
	enh.bs.ep.crt$is_tg = rep(0, length(enh.bs.ep.crt))
	
	# find overlaps with dgf
	hit.p.bg <- findOverlaps(dgf.bs,prom.bs.ep.crt,ignore.strand=TRUE,type="within")
	hit.e.bg <- findOverlaps(dgf.bs,enh.bs.ep.crt,ignore.strand=TRUE,type="within")
	
	# how many DGFs fall within the enhancers and promoters
	length(dgf.bs.reduced.bg.p <- dgf.bs[sort(queryHits(hit.p.bg))])
	length(dgf.bs.reduced.bg.e <- dgf.bs[sort(queryHits(hit.e.bg))])

	fimo.gr.bg.p <- runFimo(dgf.bs.reduced.bg.p, genome=genome, type=data_type)
	fimo.gr.bg.e <- runFimo(dgf.bs.reduced.bg.e, genome=genome, type=data_type)
	
	# Analyze FIMO results
	hoc_names = unique(c(fimo.gr.bg.p$TF,fimo.gr.p$TF))
	tf.names = rownames(df.annot.hoc[hoc_names,])
	tf.widths = df.annot.hoc[hoc_names,]$Length
	names(tf.widths) = tf.names
	res.p <- findEnrichedTFsHGtest(fimo.gr.p,fimo.gr.bg.p,dgf.bs.reduced.p,dgf.bs.reduced.bg.p,tf.names, tf.widths)
	res=res.p
	# HyperGeometric (HG) p values, one for each TF
	res.p.pv <- res$pVal
	# Enrichment factor - proportion of TF target hits divided by proportion of TF background hits
	res.p.eF <- (res$tg.hits/res$tg.lens)/((res$crt.hits+res$tg.hits)/(res$crt.lens+res$tg.lens)) 
	names(res.p.pv) = df.annot.hoc[hoc_names,]$TF
	names(res.p.eF) = df.annot.hoc[hoc_names,]$TF
	
	res.p.pv <- sort(res.p.pv)
	res.p.pv <- res.p.pv[res.p.pv!=1]
	res.p.eF <- sort(res.p.eF,decreasing=T)
	
	hoc_names = unique(c(fimo.gr.bg.e$TF,fimo.gr.e$TF))
	tf.names = rownames(df.annot.hoc[hoc_names,])
	tf.widths = df.annot.hoc[hoc_names,]$Length
	names(tf.widths) = tf.names
	res.e <- findEnrichedTFsHGtest(fimo.gr.e,fimo.gr.bg.e,dgf.bs.reduced.e,dgf.bs.reduced.bg.e,tf.names, tf.widths)
	res=res.e
	# HyperGeometric (HG) p values, one for each TF
	res.e.pv <- res$pVal
	# Enrichment factor - proportion of TF target hits divided by proportion of TF background hits
	res.e.eF <- (res$tg.hits/res$tg.lens)/((res$crt.hits+res$tg.hits)/(res$crt.lens+res$tg.lens))
	names(res.e.pv) = df.annot.hoc[hoc_names,]$TF
	names(res.e.eF) = df.annot.hoc[hoc_names,]$TF
	
	res.e.pv <- sort(res.e.pv)
	res.e.pv <- res.e.pv[res.e.pv!=1]
	res.e.eF <- sort(res.e.eF,decreasing=T)
	
	list.res[[s.cell]] = list(res.p,res.e)
	list.tf.pVal[[s.cell]] = list(res.p.pv,res.e.pv)
	list.tf.eF[[s.cell]] = list(res.p.eF,res.e.eF)
	list.tf.count[[s.cell]] = list(p.bg=sort(table(fimo.gr.bg.p$TF)),e.bg=sort(table(fimo.gr.bg.e$TF)),p.tg=sort(table(fimo.gr.p$TF)),e.tg=sort(table(fimo.gr.e$TF)))
	
	# in case the running is interrupted, we only need to re-start at last.pos, and combine the files after all going over all cells
	if(i%%20==0){
		chunk = i/20
		saveRDS(list.res[last.pos:i], file=paste0(tmp_directory,data_type,'.fimo.res.',ratio,'.',chunk,'.rds'))
		saveRDS(list.tf.pVal[last.pos:i], file=paste0(tmp_directory,data_type,'.fimo.pVal.',ratio,'.',chunk,'.rds'))
		saveRDS(list.tf.count[last.pos:i], file=paste0(tmp_directory,data_type,'.fimo.tf.count.',ratio,'.',chunk,'.rds'))
		saveRDS(list.tf.eF[last.pos:i], file=paste0(tmp_directory,data_type,'.fimo.tf.eF.',ratio,'.',chunk,'.rds'))
		last.pos = i+1
	}
	
}

## In case you need to combine the chunks - run this code (change FALSE to TRUE)
if(FALSE){
	list.res = <- list()
	list.tf.count = <- list()
	list.tf.eF = <- list()
	list.tf.pVal = <- list()
	files <- list.files(path = tmp_directory, pattern = paste0(data_type,'.fimo.'))
	 
	for(chunk in 1:(length(files)/4)){
		list.res <- c(list.res,readRDS(file=paste0(tmp_directory,data_type,'.fimo.res.',ratio,'.',chunk,'.rds')))
		list.tf.pVal <- c(list.tf.pVal,readRDS(file=paste0(tmp_directory,data_type,'.fimo.pVal.',ratio,'.',chunk,'.rds')))
		list.tf.count <- c(list.tf.count,readRDS(file=paste0(tmp_directory,data_type,'.fimo.tf.count.',ratio,'.',chunk,'.rds')))
		list.tf.eF <- c(list.tf.eF,readRDS(file=paste0(tmp_directory,data_type,'.fimo.tf.eF.',ratio,'.',chunk,'.rds')))
	}
}

saveRDS(list.res, file=paste0(tmp_directory,data_type,'.fimo.res.',ratio,'.rds'))
saveRDS(list.tf.pVal, file=paste0(tmp_directory,data_type,'.fimo.pVal.',ratio,'.rds'))
saveRDS(list.tf.count, file=paste0(tmp_directory,data_type,'.fimo.tf.count.',ratio,'.rds'))
saveRDS(list.tf.eF, file=paste0(tmp_directory,data_type,'.fimo.tf.eF.',ratio,'.rds'))

########################################################## Analyze motif finding restuls ###############################################################

list.res=readRDS(file=paste0(tmp_directory,data_type,'.fimo.res.',ratio,'.rds'))
list.tf.count=readRDS(file=paste0(tmp_directory,data_type,'.fimo.tf.count.',ratio,'.rds'))
list.tf.eF=readRDS(file=paste0(tmp_directory,data_type,'.fimo.tf.eF.',ratio,'.rds'))
list.tf.pVal = readRDS(file=paste0(tmp_directory,data_type,'.fimo.pVal.',ratio,'.rds'))

## sanity check that we actually get different number of TF hits between the target and control sets
#sapply(list.res,function(x) length(x[[1]]$crt.hits)-sum(x[[1]]$crt.hits==0))
#sapply(list.res,function(x) length(x[[1]]$tg.hits)-sum(x[[1]]$tg.hits==0))
#sapply(list.res,function(x) length(x[[2]]$crt.hits)-sum(x[[2]]$crt.hits==0))
#sapply(list.res,function(x) length(x[[2]]$tg.hits)-sum(x[[2]]$tg.hits==0))

cells = names(list.tf.count)

mat.pv.tfs.p <- matrix(0,nrow=dim(df.annot.hoc)[1],ncol=length(cells))
rownames(mat.pv.tfs.p) <- df.annot.hoc$TF
colnames(mat.pv.tfs.p) <- cells

mat.pv.tfs.e <- matrix(0,nrow=dim(df.annot.hoc)[1],ncol=length(cells))
rownames(mat.pv.tfs.e) <- df.annot.hoc$TF
colnames(mat.pv.tfs.e) <- cells

mat.eF.tfs.e <- matrix(1,nrow=dim(df.annot.hoc)[1],ncol=length(cells))
rownames(mat.eF.tfs.e) <- df.annot.hoc$TF
colnames(mat.eF.tfs.e) <- cells

mat.eF.tfs.p <- matrix(1,nrow=dim(df.annot.hoc)[1],ncol=length(cells))
rownames(mat.eF.tfs.p) <- df.annot.hoc$TF
colnames(mat.eF.tfs.p) <- cells

for ( cell in cells){
	e.tg = list.tf.count[[cell]]$e.tg
	e.tg = e.tg/sum(e.tg)
	p.tg = list.tf.count[[cell]]$p.tg
	p.tg = p.tg/sum(p.tg)
	e.eF = list.tf.eF[[cell]][[2]]
	p.eF = list.tf.eF[[cell]][[1]]
	res.p = list.tf.pVal[[cell]][[1]]
	res.e = list.tf.pVal[[cell]][[2]]
	res.p = sort(res.p)
	res.e = sort(res.e)
	gene.sym.p = names(res.p)
	gene.sym.e = names(res.e)
	print(paste0(cell, ':'))
	print("Promoter TFs:")
	print(res.p)
	mat.pv.tfs.p[names(res.p),cell] = (-1)*log10(res.p)	
	print("Enhancer TFs:")
	print(res.e)
	mat.pv.tfs.e[names(res.e),cell] = (-1)*log10(res.e)
	mat.eF.tfs.e[names(res.e),cell] = e.eF[names(res.e)]
	mat.eF.tfs.p[names(res.p),cell] = p.eF[names(res.p)]
}

saveRDS(mat.pv.tfs.e, file=paste0(tmp_directory,data_type,'.mat.pv.tfs.e.rds'))
saveRDS(mat.pv.tfs.p, file=paste0(tmp_directory,data_type,'.mat.pv.tfs.p.rds'))
saveRDS(mat.eF.tfs.e, file=paste0(tmp_directory,data_type,'.mat.eF.tfs.e.rds'))
saveRDS(mat.eF.tfs.p, file=paste0(tmp_directory,data_type,'.mat.eF.tfs.p.rds'))

## Specificity score of the TFs' enrichment factors

## Specificity source functions
source(paste0(script_directory,'FUNCTIONS_specificity_score.R'))

Sc.e = rep(list(NULL),length(cells))
Sc.p = rep(list(NULL),length(cells))
ranks.p = rep(length(cells),length(cells))
ranks.e = rep(length(cells),length(cells))

names(Sc.e) <- names(Sc.p) <- cells
names(ranks.e) <- names(ranks.p) <- cells

num.tf.e = rep(0,length(cells))
num.tf.p = rep(0,length(cells))
names(num.tf.e) <- names(num.tf.p) <- cells

for ( cell in cells){
	flag=TRUE
	try({
	x <- mat.eF.tfs.e
	y <- mat.pv.tfs.e
	sub_tfs <- names(which(y[,cell]>(-log10(pv.th))))
	x <- as.matrix(x[sub_tfs,,drop=F])
	num.tf.e[cell] = dim(x)[1]
	Sc.e[[cell]] = unlist(mclapply(cells,function(y) cellSpecificityScore(y,x),mc.cores=10))
	names(Sc.e[[cell]]) = cells
	ranks.e[cell] = which(names(sort(Sc.e[[cell]],decreasing=T))==cell)
	
	x <- mat.eF.tfs.p
	y <- mat.pv.tfs.p
	sub_tfs <- names(which(y[,cell]>(-log10(pv.th))))
	x <- as.matrix(x[sub_tfs,,drop=F])
	num.tf.p[cell] = dim(x)[1]
	Sc.p[[cell]] = unlist(mclapply(cells,function(y) cellSpecificityScore(y,x),mc.cores=10))
	names(Sc.p[[cell]]) = cells
	ranks.p[cell] = which(names(sort(Sc.p[[cell]],decreasing=T))==cell)
	flag=FALSE
	})
	if(flag)
		print(cell)
}

saveRDS(Sc.e, file=paste0(tmp_directory,data_type,'.TF.specificity.enhancers.rds'))
saveRDS(Sc.p, file=paste0(tmp_directory,data_type,'.TF.specificity.promoters.rds'))


