##################################################################################################################################
################################# First step modeling (part 2) - Finding robust promoter models ##################################
##################################################################################################################################

########################################## Note: Only works via R unix platform ##################################################

## This is a line by line running tutorial. In order to to run a subset of lines please select them using shift and run using ctrl+R

## Please follow the instruction given in http://acgt.cs.tau.ac.il/ct-focs/tutorial.html

# load libraries - make sure you pre-installed these packages
cran_libs = c('pscl','MASS','parallel','AUC','glmnet','RColorBrewer','ggplot2')
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
# data types - fantom, encode
data_type <- 'encode'

#load code
source(paste(script_directory,'FUNCTIONS_reg_functions.R',sep=''))

# Enhancer genomic positions
enh.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.enh.pos.rds',sep=''))
# Promoters positions
prom.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.prom.pos.rds',sep=''))
# Enhancer count matrix:
Me = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.enh.rle.rds',sep='')))
# Promoter count matrix:
Mg = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.prom.rle.rds',sep='')))

# Sample annotations
sample_annot <- readRDS(paste(data_directory,data_type,'.sample.annot.rds',sep=''))
# Promoter to closest enhnacers in fixed window:
g_to_e = readRDS(paste(data_directory,data_type,'/',data_type,'.cand.enh.rds',sep=''))


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


auroc_m <- as.numeric(readRDS(file=paste(tmp_directory,data_type,'.auroc_m.rds',sep='')))
corr_m <- as.numeric(readRDS(file=paste(tmp_directory,data_type,'.corr_m.rds',sep='')))
corr_pvals_m <- as.numeric(readRDS(file=paste(tmp_directory,data_type,'.corr_pvals_m.rds',sep='')))
preds_m <- readRDS(file=paste(tmp_directory,data_type,'.preds_m.rds',sep=''))
auroc_pvals_m <- readRDS(file=paste(tmp_directory,data_type,'.auroc_pvals_m.rds',sep=''))
r_square_lm <- readRDS(file=paste(tmp_directory,data_type,'.r_square_lm.rds',sep=''))
r2<-r_square_lm

gene2num_positive_samples = c()
if(data_type=='encode')
	gene2num_positive_samples = apply(Mg_fpkm>1,1,sum,na.rm=T)
if(data_type=='fantom5')
	gene2num_positive_samples = apply(Mg_fpkm>0,1,sum,na.rm=T)

names(gene2num_positive_samples) = rownames(Mg)

# Look at the OLS results: obtain estimation for the FDR and power of using R.squared values
binary_validation_qvals = p.adjust(auroc_pvals_m,method=method)
binary_validation_qvals[is.na(binary_validation_qvals)] = 1
binary_validation_qvals = set_min_greater_than_zero(binary_validation_qvals)
level_validation_qvals = p.adjust(corr_pvals_m,method=method)
level_validation_qvals[is.na(level_validation_qvals)] = 1
level_validation_qvals = set_min_greater_than_zero(level_validation_qvals)

# Divide promoter models into groups based on the FDR threshold
gene_val_groups  = list()
# Promoter models that passed both 'binary' and 'expression level' validations
gene_val_groups[["Both"]] = names(which(binary_validation_qvals<=fdr_th & level_validation_qvals<=fdr_th))
# Promoter models that passed only the binary validation
gene_val_groups[["Binary only"]] = names(which(binary_validation_qvals<=fdr_th & !level_validation_qvals<=fdr_th))
# Promoter models that passed only the expression level validation
gene_val_groups[["Level only"]] = names(which(!binary_validation_qvals<=fdr_th & level_validation_qvals<=fdr_th))
# Promoter models that did not pass any validation
gene_val_groups[["None"]] = names(which(!binary_validation_qvals<=fdr_th & !level_validation_qvals<=fdr_th))

saveRDS(gene_val_groups,file=paste0(tmp_directory,data_type,'.gene_val_groups.rds'))

################################################ Performance plots ##########################################################

# Figure 1.A-B
# Pie chart of the number of models passed each group (Both, Binary only, Level only)
# Boxplot of the number of positive samples each model contains across the groups
# You may want to change some parameters, e.g. ylim/width/height, to adjust the plots

png(paste(tmp_directory,data_type,'.pie_box.png',sep=''),width=1000,height=1200)
par(mfrow=c(2,1),mar=c(6,6,2,1),cex=2.0)
marker = list(color = brewer.pal(4, "Pastel1"))
#cols = topo.colors(4)
cols = marker$color
(group_sizes <- sapply(gene_val_groups,length))
names(cols) = names(group_sizes)
labls = paste("\n\n",names(group_sizes),group_sizes,"\n",sep="\n")
labls[1] = paste(names(group_sizes)[1],group_sizes[1],"\n",sep="\n")
labls[2] = paste("       \n\n",names(group_sizes)[2],group_sizes[2],"\n    ",sep="\n")
pie(group_sizes,labels=labls,col=cols,lwd=2,cex=1,radius=0.9)
expr_by_group = lapply(gene_val_groups,function(x,y)y[x],y=gene2num_positive_samples)
boxplot(expr_by_group,ylab = paste0("Number of positive samples","\n","[RLE>1]"),las=2,col=cols,cex.lab=1.2,cex.names=1.2)

dev.off()


# Figure 1.C
# R.squared without CV vs. R.squared with CV
# You may want to change some parameters, e.g. ylim/width/height, to adjust the plo

r_cv <- sapply(seq_len(dim(Mg_fpkm)[1]),function(i) getRsquared(i,Mg_fpkm[i,],preds_m[rownames(Mg_fpkm)[i],]))

png(paste(tmp_directory,data_type,'.r_cv_vs_r2.png',sep=''),width=1000,height=1000)
par(mar=c(6,6,1,1),cex=1.8)
cols = rep("black",length(r2))
cols[r_square_lm>=0.5] = "red"
cols[r_square_lm>=0.5 & r_cv >= 0.25] = "blue"
plot(r_square_lm,r_cv,xlab = expression("R"^2), ylab = expression("R"[CV]^2), cex.lab=1.6,cex.axis=1.6,col=cols,mex=1.15)
labels <- c(expression(paste("R"^2>=0.5," & ","R"[CV]^2<0.25)),expression(paste("R"^2>=0.5," & ","R"[CV]^2>=0.25)))
legend("topleft",legend=labels,col=c("red","blue"),pch=21,cex=1.4)
table(cols)
sum(cols=="red")/sum(cols!="black")
dev.off()



