##################################################################################################################################
################################# First step modeling (part 1) - Finding robust promoter models ##################################
##################################################################################################################################


########################################## Note: Only works via R unix platform ##################################################

## This is a line by line running tutorial. In order to to run a subset of lines please select them using shift and run using ctrl+R

## Please follow the instruction given in http://acgt.cs.tau.ac.il/ct-focs/tutorial.html

# load libraries - make sure you pre-installed these packages
cran_libs = c('pscl','MASS','parallel','AUC','RColorBrewer','ggplot2')
bioconductor_libs = c('GenomicRanges')
for (p in (c(bioconductor_libs,cran_libs))){library(p,character.only=T)}


# Variables
k <- 10 #k closest enhancers to each gene
win.size <- 5*10**5 #maximum window size upstream/downstream to select genes with at least 10 enhancers
# Number of parallel cores to run on - use detectCores(all.tests = FALSE, logical = TRUE) to identify the number of cores
mc.cores <- 40
# FDR correction method
method='BY'
# FDR threshold
fdr_th = 0.1
# Regression type (1:ZINB, 2:OLS, 3:GLM.NB)
reg_type <- 2

# load objects and read files
data_directory = 'data/'
script_directory = 'scripts/'
tmp_directory = 'tmp/'
# data types - fantom5, encode
data_type <- 'encode'


# Enhancer genomic positions
enh.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.enh.pos.rds',sep=''))
# Promoters positions
prom.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.prom.pos.rds',sep=''))
# Enhancer RLE matrix:
Me = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.enh.rle.rds',sep='')))
# Promoter RLE matrix:
Mg = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.prom.rle.rds',sep='')))

# Sample annotations
sample_annot <- readRDS(paste(data_directory,data_type,'/',data_type,'.sample.annot.rds',sep=''))
# Promoter to closest enhnacers in fixed window:
g_to_e = readRDS(paste(data_directory,data_type,'/',data_type,'.cand.enh.rds',sep=''))

#load code
source(paste(script_directory,'FUNCTIONS_reg_functions.R',sep=''))

tmpMe <- Me
tmpMg <- Mg[names(g_to_e),]
col2type <- sample_annot[,'col2type',drop=F]

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

pos.th <- 1 # Promoter activity threshold of positive samples
num.pos.th <- 3 # how many samples should have at least pos.th RLE promoter activity

if(data_type=='fantom5')
	pos.th <- 0

col2type = as.character(col2type[,1])
# list of arguments
arg_list <- list(Mg,Mg_normalized,Me,Me_normalized,g_to_e,offs,Mg_fpkm,col2type)
res <- mclapply(seq_len(nrow(Mg)),function(j) apply_cv_test_by_single_type(j,arg_list,k,func = run_lm_test,pos.th=pos.th,num.pos.th=num.pos.th),mc.cores = mc.cores)

auroc_m <- do.call(c,sapply(res,function(x) x[1]))
corr_m <- do.call(c,sapply(res,function(x) x[2]))
corr_pvals_m <- do.call(c,sapply(res,function(x) x[3]))
preds_m <- do.call(rbind,sapply(res,function(x) x[4][[1]]))

names(auroc_m) = rownames(Mg)
names(corr_m) = rownames(Mg)
names(corr_pvals_m) = rownames(Mg)
rownames(preds_m) = rownames(Mg)

saveRDS(auroc_m,file=paste(tmp_directory,data_type,'.auroc_m.rds',sep=''))
saveRDS(corr_m,file=paste(tmp_directory,data_type,'.corr_m.rds',sep=''))
saveRDS(corr_pvals_m,file=paste(tmp_directory,data_type,'.corr_pvals_m.rds',sep=''))
saveRDS(preds_m,file=paste(tmp_directory,data_type,'.preds_m.rds',sep=''))


# Get ROC scores p-values
auroc_pvals_m = c()
for(i in 1:nrow(Mg)){
  if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",i))}
  y_fpkm = Mg[i,]
  if(length(g_to_e[[i]]) < k || any(is.na(y_fpkm)) || sum(y_fpkm > pos.th) < num.pos.th){
    res = c(NA)
    auroc_pvals_m = c(auroc_pvals_m,res)
    next
  }
  if(sum(y_fpkm <= pos.th) ==0){
    res = c(NA)
    auroc_pvals_m = c(auroc_pvals_m,res)
    next
  }
  pos_inds = y_fpkm > pos.th
  neg_inds = !pos_inds
  curr_pvals = c()
  curr_preds = as.numeric(preds_m[i,])
  x1 = curr_preds[pos_inds]
  x2 = curr_preds[neg_inds]
  na_ids <- is.na(x1)
  if(sum(na_ids)!=0) x1 <- x1[!na_ids]
  na_ids <- is.na(x2)
  if(sum(na_ids)!=0) x2 <- x2[!na_ids]	
  if(length(x1)==0 || length(x2)==0){
    auroc_pvals_m = c(auroc_pvals_m,NA)
    next
  }
  curr_pvals = wilcox.test(x1,x2)$p.value
  auroc_pvals_m = c(auroc_pvals_m,curr_pvals)			
}
names(auroc_pvals_m) = rownames(Mg)
saveRDS(auroc_pvals_m,file=paste(tmp_directory,data_type,'.auroc_pvals_m.rds',sep=''))


r_square_lm <- c()
for(j in 1:nrow(Mg)){
  if(j%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",j))}
  p_name <- rownames(Mg)[j]
  y_fpkm <- Mg_fpkm[p_name,]
  if(length(g_to_e[[p_name]]) < k || any(is.na(y_fpkm)) || sum(y_fpkm > pos.th) < num.pos.th){
    r_square_lm <- c(r_square_lm,NA)
    next
  }
  regr_data_normalized = fetch_regr_data(p_name,Mg,Me,g_to_e)
  x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
  #x_n <- x_n[rownames(primary_cells),]; y_n <- y_n[rownames(primary_cells)]
  if(ncol(x_n)<k){
    r_square_lm <- c(r_square_lm,NA)
    next
  }
  if(ncol(x_n)>k) x_n <- x_n[,1:k]
  res <- run_lm_test(x_n,y_n)
  r_square_lm <- c(r_square_lm,res$r.squared)
  
}
names(r_square_lm) <- rownames(Mg)
saveRDS(r_square_lm,file=paste(tmp_directory,data_type,'.r_square_lm.rds',sep=''))

