##################################################################################################################################
#################################################### Paper's Figures #############################################################
##################################################################################################################################

########################################## Note: Only works via R unix platform ##################################################

## This is a line by line running tutorial. In order to to run a subset of lines please select them using shift and run using ctrl+R

## Please follow the instruction given in http://acgt.cs.tau.ac.il/ct-focs/tutorial.html

# load libraries - make sure you pre-installed these packages
cran_libs = c('S4Vectors','RColorBrewer','ggplot2','gplots')
bioconductor_libs = c('GenomicRanges','biomaRt')
for (p in (c(bioconductor_libs,cran_libs))){library(p,character.only=T)}

# load objects and read files
data_directory = 'data/'
script_directory = 'scripts/'
tmp_directory = 'tmp/'
# data types - fantom5, encode
data_type <- 'encode'

pv.th = 0.05 # HyperGeometric p value threshold of significant enriched transcription factor (TF)

source(paste0(script_directory,'FUNCTIONS_specificity_score.R'))


s.cell <- 'GM12878'

# Enhancer genomic positions
enh.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.enh.pos.rds',sep=''))
# Promoters positions
prom.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.prom.pos.rds',sep=''))
# Enhancer RLE matrix:
Me = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.enh.rle.log2.rds',sep='')))
# Promoter RLE matrix:
Mg = as.matrix(readRDS(paste(data_directory,data_type,'/',data_type,'.prom.rle.log2.rds',sep='')))
# Sample annotations
sample_annot <- readRDS(paste(data_directory,data_type,'/',data_type,'.sample.annot.rds',sep=''))
col2type <- sample_annot[,'col2type',drop=F]

# Predicted cell type-specific EP links from CT-FOCS
ep_links = readRDS(file=paste0(data_directory,data_type,'/',data_type,'.EP_links.cell.rds'))

# load GM12878 ChIA-PET connected loop sets (CLSs)
chia_pet_CLS = readRDS(file=paste0(data_directory,data_type,'/',data_type,'.chia_pet.connected.loop.sets.',s.cell,'.rds'))

# Load ENCODE 72 cell line gene expression data from Sheffield et al., Genome Research (2013)
exp.data <- readRDS(file=paste0(data_directory,"Sheffield_GR_2013.exp.data.rds"))
exp.data <- log2(exp.data+1)

e.active.cell <- ep_links[[s.cell]]
p.active.cell <- names(e.active.cell)
length(p.active.cell)

tmp <- e.active.cell
tmp[sapply(tmp,is.null)] <- NULL

length(unlist(tmp)) # num of EP links in s.cell
length(unique(unlist(tmp))) # num of linked enhancers in s.cell
length(tmp) # num of linked promoters in s.cell

ep.list = tmp

res <- createTPepSet(ep.list,chia_pet_CLS)
true_ep <- res[[1]]

num_sup_per_p = sapply(true_ep,sum)
num_sup_per_p[which.max(num_sup_per_p)]
prom.bs[names(num_sup_per_p)[which(num_sup_per_p>=3)]]$sym_id

x <- createDijMatLog(ep.list,Me,Mg,col2type)
colnames(x) <- col2type[colnames(Mg),'col2type']
y <- x

## Fig. 3A-B
#initiate cols with all black
cols <- rep('red', nrow(x))
cols[unlist(true_ep,use.names=F)] <- 'yellow'

colsCol <- rep('black', ncol(x))
colsCol[grepl('^GM', colnames(x))] <- 'blue'
colsCol[grepl('^Th', colnames(x))] <- 'green'
colsCol[which(colnames(x)=='GM04504')]='black' # not a lymphocyte
colsCol[which(colnames(x)=='GM04503')]='black' # not a lymphocyte
colsCol[which(colnames(x)==s.cell)] <- 'red'

cols = factor(cols,levels=c("yellow","red"))
ix = sort.int(cols,index.return=TRUE)$ix
x <- x[ix,]
cols <- as.character(cols[ix])

# Calculate cell type specificity score based on the EP signals
Sc = mclapply(unique(colnames(y)), function(cell) cellSpecificityScore(cell,y), mc.cores=10)
names(Sc) = unique(colnames(y))
Sc = unlist(Sc)
Sc = sort(Sc,decreasing=T)
Sc.col = rep("grey",length(Sc))
names(Sc.col) <- names(Sc)
Sc.col[grepl('^GM', names(Sc.col))] <- 'blue'
Sc.col[which(names(Sc.col)=='GM04504')]='grey' # not a lymphocyte
Sc.col[which(names(Sc.col)=='GM04503')]='grey' # not a lymphocyte
Sc.col[grepl('^Th', names(Sc.col))] <- 'green'
Sc.col[s.cell] = "red"

Sc.cex = rep(1.5,length(Sc))
names(Sc.cex) <- names(Sc)
Sc.cex[s.cell] = 2.5

Sc.ep = Sc; Sc.cex.ep = Sc.cex; Sc.col.ep = Sc.col

pdf(file=paste0("pics/",data_type,'.',s.cell,".heatmap.ep.pdf"))
test <- heatmap.2(x,
                  trace="none",
                  scale="row",
                  Rowv=FALSE,
                  labRow = NA,
                  labCol = NA,
                  distfun = dist,
                  symbreaks = TRUE,
                  colRow = cols,
                  RowSideColors = cols,
                  colCol = colsCol,
                  key.title = NA,
                  ColSideColors = colsCol,
                  col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
                  margins=c(4,7),breaks = seq(-3, 3, length.out = 101)
)
dev.off()

## Gene Expression (GE) heatmap

length(gene.sym <- unique(unlist(prom.bs[names(ep.list)]$sym_id)))
exp.data.tmp <- exp.data[which(!is.na(match(rownames(exp.data),gene.sym))),]
x = as.matrix(exp.data.tmp)
y = x

colsCol <- rep('black', ncol(x))
colsCol[grepl('^GM', colnames(x))] <- 'blue'
colsCol[which(colnames(x)=='Th1')]='green'
colsCol[which(colnames(x)=='Th2')]='green'
colsCol[which(colnames(x)=='CLL')]='green'
colsCol[which(colnames(x)=='GM04504')]='black' # not a lymphocyte
colsCol[which(colnames(x)=='GM04503')]='black' # not a lymphocyte
colsCol[which(colnames(x)==s.cell)] <- 'red'


# Calculate cell type specificity score based on the GE
Sc = mclapply(unique(colnames(y)), function(cell) cellSpecificityScore(cell,y), mc.cores=10)
names(Sc) = unique(colnames(y))
Sc = unlist(Sc)
Sc = sort(Sc,decreasing=T)
Sc.col = rep("grey",length(Sc))
names(Sc.col) <- names(Sc)
Sc.col[grepl('^GM', names(Sc.col))] <- 'blue'
Sc.col[which(names(Sc.col)=='GM04504')]='grey' # not a lymphocyte
Sc.col[which(names(Sc.col)=='GM04503')]='grey' # not a lymphocyte
Sc.col[grepl('^Th', names(Sc.col))] <- 'green'
Sc.col[which(names(Sc.col)=='CLL')]='green'
Sc.col[s.cell] = "red"

Sc.cex = rep(1.5,length(Sc))
names(Sc.cex) <- names(Sc)
Sc.cex[s.cell] = 2.5

Sc.ge = Sc; Sc.cex.ge = Sc.cex; Sc.col.ge = Sc.col


pdf(file=paste0("pics/",data_type,'.',s.cell,".heatmap.ge.pdf"))
test <- heatmap.2(x,
                  trace="none",
                  scale="row",
                  Rowv=FALSE,
                  #Colv=as.dendrogram(hc),
                  labRow = NA,
                  labCol = NA,
                  key.title = NA,
                  distfun = dist,
                  symbreaks = TRUE,
                  colCol = colsCol,
                  ColSideColors = colsCol,
                  col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
                  margins=c(6,7),breaks = seq(-3, 3, length.out = 101)
)

dev.off()

## Fig. 3C-D


pdf(file=paste0("pics/",data_type,'.',s.cell,".ep.spec.pdf"))
plot(1:length(Sc.ep),range01(Sc.ep),type='p',pch=19,col=Sc.col.ep,ylab="Cell type specificity score",xlab="Cell types",cex=Sc.cex.ep)
legend('topright',legend=c(s.cell,"Other B-cells","T-cells","Other cell types"),col=c("red","blue","green","grey"),pch=19)
dev.off()

pdf(file=paste0("pics/",data_type,'.',s.cell,".GE.spec.pdf"))
plot(1:length(Sc.ge),range01(Sc.ge),type='p',pch=19,col=Sc.col.ge,ylab="Cell type specificity score",xlab="Cell types",cex=Sc.cex.ge)
legend('topright',legend=c(s.cell,"Other B-cells","T-cells","Other cell types"),col=c("red","blue","green","grey"),pch=19)
dev.off()

## Fig. 4A

data_type = 'fantom5'

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
col2type <- sample_annot[,'col2type',drop=F]

# Predicted cell type-specific EP links from CT-FOCS
ep_links = readRDS(file=paste0(data_directory,data_type,'/',data_type,'.EP_links.cell.rds'))

## CT-FOCS

e.active.cell <- ep_links[[tolower(s.cell)]]
p.active.cell <- names(e.active.cell)
length(p.active.cell)

tmp <- e.active.cell
tmp[sapply(tmp,is.null)] <- NULL

length(unlist(tmp)) # num of EP links in s.cell
length(unique(unlist(tmp))) # num of linked enhancers in s.cell
length(tmp) # num of linked promoters in s.cell

ep.list = tmp

length(gene.sym <- unique(unlist(prom.bs[names(ep.list)]$sym_id)))
exp.data.tmp <- exp.data[which(!is.na(match(rownames(exp.data),gene.sym))),]
x = as.matrix(exp.data.tmp)
y = x


colsCol <- rep('black', ncol(x))
colsCol[grepl('^GM', colnames(x))] <- 'blue'
colsCol[which(colnames(x)=='Th1')]='green'
colsCol[which(colnames(x)=='Th2')]='green'
colsCol[which(colnames(x)=='CLL')]='green'
colsCol[which(colnames(x)=='GM04504')]='black' # not a lymphocyte
colsCol[which(colnames(x)=='GM04503')]='black' # not a lymphocyte
colsCol[which(colnames(x)==s.cell)] <- 'red'


pdf(file=paste0("pics/",data_type,'.',s.cell,".heatmap.ge.ct_focs.pdf"))
test <- heatmap.2(x,
                  trace="none",
                  scale="row",
                  Rowv=FALSE,
                  #Colv=as.dendrogram(hc),
                  labRow = NA,
                  labCol = NA,
                  key.title = NA,
                  distfun = dist,
                  symbreaks = TRUE,
                  colCol = colsCol,
                  ColSideColors = colsCol,
                  col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
                  margins=c(6,7),breaks = seq(-3, 3, length.out = 101)
)

dev.off()

## JEME
score = 0.3
prom.gr.jeme = readRDS(file=paste0(data_directory,data_type,'.prom.pos.jeme.rds'))
enh.gr.jeme = readRDS(file=paste0(data_directory,data_type,'.enh.pos.jeme.rds'))

enh.gr.jeme.cell <- enh.gr.jeme[unlist(mclapply(tolower(enh.gr.jeme$tissue),function(x) (tolower(s.cell) %in% x),mc.cores=10))]
prom.gr.jeme.cell <- prom.gr.jeme[unlist(mclapply(tolower(prom.gr.jeme$tissue),function(x) (tolower(s.cell) %in% x),mc.cores=10))]

match.id <- which(prom.gr.jeme.cell$score>=score)
enh.gr.jeme.cell.tmp <- enh.gr.jeme.cell[match.id]
prom.gr.jeme.cell.tmp <- prom.gr.jeme.cell[match.id]
ep.int <- paste0(names(enh.gr.jeme.cell.tmp),'##',names(prom.gr.jeme.cell.tmp))
ep.int <- unique(ep.int)
e.name <- unlist(strsplit(ep.int,"##"))[seq(1,2*length(ep.int),2)]
p.name <- unlist(strsplit(ep.int,"##"))[seq(2,2*length(ep.int),2)]
enh.gr.jeme.cell.tmp <- enh.gr.jeme.cell.tmp[e.name]
prom.gr.jeme.cell.tmp <- prom.gr.jeme.cell.tmp[p.name]
  
length(enh.gr.jeme.cell.tmp) #num of EP
length(unique(enh.gr.jeme.cell.tmp)) #num of E
length(unique(prom.gr.jeme.cell.tmp)) #num of P

length(gene.sym <- unique(unlist(prom.gr.jeme.cell.tmp$sym_id)))
exp.data.tmp <- exp.data[which(!is.na(match(rownames(exp.data),gene.sym))),]
x = as.matrix(exp.data.tmp)
y = x


colsCol <- rep('black', ncol(x))
colsCol[grepl('^GM', colnames(x))] <- 'blue'
colsCol[which(colnames(x)=='Th1')]='green'
colsCol[which(colnames(x)=='Th2')]='green'
colsCol[which(colnames(x)=='CLL')]='green'
colsCol[which(colnames(x)=='GM04504')]='black' # not a lymphocyte
colsCol[which(colnames(x)=='GM04503')]='black' # not a lymphocyte
colsCol[which(colnames(x)==s.cell)] <- 'red'

## Warning - the matrix is of size 4,048 x 112 - this will take a long time to create the heatmap

pdf(file=paste0("pics/",data_type,'.',s.cell,".heatmap.ge.jeme.pdf"))
test <- heatmap.2(x,
                  trace="none",
                  scale="row",
                  Rowv=FALSE,
                  #Colv=as.dendrogram(hc),
                  labRow = NA,
                  labCol = NA,
                  key.title = NA,
                  distfun = dist,
                  symbreaks = TRUE,
                  colCol = colsCol,
                  ColSideColors = colsCol,
                  col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
                  margins=c(6,7),breaks = seq(-3, 3, length.out = 101)
)

dev.off()

## Fig. 4B

Sc_exp_ranks_ct_focs=readRDS(file=paste0(tmp_directory,"fantom.ge.sc.scores.ct_focs.rds"))
Sc_exp_ranks_jeme=readRDS(file=paste0(tmp_directory,"fantom.ge.sc.scores.jeme.rds"))

cells = intersect(names(Sc_exp_ranks_jeme),names(Sc_exp_ranks_ct_focs))

# The results are shown as table for 4 cel types: GM12878, HepG2, K562, and MCF-7
df_exp = data.frame("CT_FOCS"=Sc_exp_ranks_ct_focs[cells],'JEME'=Sc_exp_ranks_jeme[cells])
df_exp[order(df_exp$CT_FOCS),]

## Fig. 4C

Sc_ranks_jeme=readRDS(file=paste0(tmp_directory,"fantom.ep.sc.scores.jeme.rds"))
Sc_ranks_ct_focs=readRDS(file=paste0(tmp_directory,"fantom.ep.sc.scores.ct_focs.rds"))

cells = intersect(names(Sc_ranks_jeme),names(Sc_ranks_ct_focs))

df_ep = data.frame("CT_FOCS"=Sc_ranks_ct_focs[cells],'JEME'=Sc_ranks_jeme[cells])
head(df_ep[order(df_ep$CT_FOCS),],20)

df = data.frame('Sc_rank'=c(Sc_ranks_ct_focs[cells],Sc_ranks_jeme[cells]),'Method'=c(rep('CT-FOCS',length(cells)),rep('JEME',length(cells))))
cdf = data.frame(Method=c('CT-FOCS','JEME'), Sc_mean=c(mean(Sc_ranks_ct_focs[cells]),mean(Sc_ranks_jeme[cells])))

pdf(paste0("pics/",data_type,".cellType.spec.score.ranks.ct_focs_jeme.density.pdf"), width=10, height=5)
p1 = ggplot(data=df,aes(x=Sc_rank,fill=Method, colour=Method))+ geom_density(alpha=0.3)+ geom_vline(data=cdf, aes(xintercept=Sc_mean,  colour=Method),
               linetype="dashed", size=1)+ xlab("Cell Type Specificity Score Rank (out of 328 cell types)") + scale_fill_manual(values=c("#1a85ff","#d41159")) + scale_colour_manual(values=c("#1a85ff","#d41159"))# + theme(legend.pos = "none")
p1 + theme(axis.text.x = element_text( color="black", 
                           size=14),
          axis.text.y = element_text( color="black", 
                           size=14))
dev.off()

# One-sided Wicoxon paired test
(p.v=wilcox.test(Sc_ranks_ct_focs[cells],Sc_ranks_jeme[cells],paired=TRUE,alternative="less")$p.value)

## Fig. 4D
s.cell = 'GM12878'

rel_specif_cell_jeme = readRDS(file=paste0(tmp_directory,'rel_specif_cell_jeme.',s.cell,'.rds'))
rel_specif_cell_ct_focs = readRDS(file=paste0(tmp_directory,'rel_specif_cell_ct_focs.',s.cell,'.rds'))

cells = intersect(names(rel_specif_cell_ct_focs),names(rel_specif_cell_jeme))

df = data.frame('Rel_Exp'=log2(c(rel_specif_cell_ct_focs[cells],rel_specif_cell_jeme[cells])),'Method'=c(rep('CT-FOCS',length(cells)),rep('JEME',length(cells))))
cdf = data.frame(Method=c('CT-FOCS','JEME'), Rel_exp_mean=c(mean(log2(rel_specif_cell_ct_focs[cells])),mean(log2(rel_specif_cell_jeme[cells]))))

pdf(paste0("pics/",data_type,'.',s.cell,".relative_precision.ctfocs.jeme.pdf"), width=10, height=5)

p2 = ggplot(data=df,aes(y=Rel_Exp,x=Method, fill=Method, colour=Method))+ geom_boxplot(alpha=0.5, outlier.shape=16, outlier.size=1, notch=TRUE) + ylab("GM12878 Enrichment [log2[fold change]]")
p2 = p2 + geom_jitter(shape=16, position=position_jitter(0.2)) + scale_fill_manual(values=c("#1a85ff","#d41159")) + scale_colour_manual(values=c("#1a85ff","#d41159"))

dev.off()

# One-sided Wicoxon paired test
(p.v=wilcox.test(rel_specif_cell_ct_focs[cells],rel_specif_cell_jeme[cells],paired=TRUE,alternative='greater')$p.value)

## Fig. 5A-B
s.cell = 'GM12878'
data_type = 'encode'

mat.pv.tfs.e <- readRDS(file=paste0(tmp_directory,data_type,'.mat.pv.tfs.e.rds'))
mat.pv.tfs.p <- readRDS(file=paste0(tmp_directory,data_type,'.mat.pv.tfs.p.rds'))
mat.eF.tfs.e <- readRDS(file=paste0(tmp_directory,data_type,'.mat.eF.tfs.e.rds'))
mat.eF.tfs.p <- readRDS(file=paste0(tmp_directory,data_type,'.mat.eF.tfs.p.rds'))

# Promoters
x <- mat.eF.tfs.p
y <- mat.pv.tfs.p
sub_tfs <- apply(y[,s.cell,drop=F],1,function(row) any(row>(-log10(pv.th))))
x <- t(as.matrix(x))
x <- t(scale(x))
x <- x[sub_tfs,,drop=F]

filename = paste0("pics/",data_type,".TFs.prom.",s.cell,".pdf")

pdf(filename,width=7, height=6, pointsize=8)
test <- heatmap.2(x,
          trace="none",
          scale="none",
          distfun = function(x) as.dist(1-cor(t(x))),
	    hclustfun = function(x) hclust(x, method="ave"),
          colCol = NA,
	    key.title = NA,
	    key.size=3.5,
	    key.xlab='Enrichment Factor (Z-score)',
 	    col = colorRampPalette(c('#fef0d9','#fdbb84','#b30000'))(20),
	    margins=c(10,7),breaks = seq(-2, 2, length.out = 21)
)

dev.off()

# Enhancers

x <- mat.eF.tfs.e
y <- mat.pv.tfs.e
sub_tfs <- apply(y[,s.cell,drop=F],1,function(row) any(row>(-log10(pv.th))))
x <- t(as.matrix(x))
x <- t(scale(x))
x <- x[sub_tfs,,drop=F]

filename = paste0("pics/",data_type,".TFs.enh.",s.cell,".pdf")


pdf(filename,width=7, height=6, pointsize=10)

test <- heatmap.2(x,
          trace="none",
          scale="none",
          distfun = function(x) as.dist(1-cor(t(x))),
	    hclustfun = function(x) hclust(x, method="ave"),
          colCol = NA,
	    key.title = NA,
	    key.size=3.5,
	    key.xlab='Enrichment Factor (Z-score)',
 	    col = colorRampPalette(c('#fef0d9','#fdbb84','#b30000'))(20),
	    margins=c(10,7),breaks = seq(-2, 2, length.out = 21)
)

dev.off()

## Fig. 5C-D

Sc.e <- readRDS(file=paste0(tmp_directory,data_type,'.TF.specificity.enhancers.rds'))
Sc.p <- readRDS(file=paste0(tmp_directory,data_type,'.TF.specificity.promoters.rds'))

# Enhancers
Sc = sort(Sc.e[[s.cell]],decreasing=T)
Sc.col = rep("grey",length(Sc))
names(Sc.col) <- names(Sc)
Sc.col[grepl('^GM', names(Sc.col))] <- 'blue'
Sc.col[grepl('^Th', names(Sc.col))] <- 'green'
Sc.col[s.cell] = "red"

Sc.cex = rep(2.5,length(Sc))
names(Sc.cex) <- names(Sc)
Sc.cex[s.cell] = 3.5
Sc.cex[tolower(s.cell)] = 3.5

pdf(file=paste0("pics/",data_type,'.',s.cell,".tf.spec.enh.pdf"))
plot(1:length(Sc),range01(Sc),type='p',pch=19,col=Sc.col,ylab="Cell type specificity score",xlab="Cell types",cex=Sc.cex, cex.axis=1.5)
legend('topright',legend=c(s.cell,"Other B-cells","T-cells","Other cell types"),col=c("red","blue","green","grey"),pch=19)
dev.off()

# Promoters
Sc = sort(Sc.p[[s.cell]],decreasing=T)
Sc.col = rep("grey",length(Sc))
names(Sc.col) <- names(Sc)
Sc.col[grepl('^GM', names(Sc.col))] <- 'blue'
Sc.col[grepl('^Th', names(Sc.col))] <- 'green'
Sc.col[s.cell] = "red"

Sc.cex = rep(2.5,length(Sc))
names(Sc.cex) <- names(Sc)
Sc.cex[s.cell] = 3.5
Sc.cex[tolower(s.cell)] = 3.5

pdf(file=paste0("pics/",data_type,'.',s.cell,".tf.spec.prom.pdf"))
plot(1:length(Sc),range01(Sc),type='p',pch=19,col=Sc.col,ylab="Cell type specificity score",xlab="Cell types",cex=Sc.cex, cex.axis=1.5)
legend('topright',legend=c(s.cell,"Other B-cells","T-cells","Other cell types"),col=c("red","blue","green","grey"),pch=19)
dev.off()

################################################ Do Not Run ############################################################

## Supplementary Figures - not all of them

### Figures

## Fig. S1A - Histogram of the number of EP links per cell type
pdf(file='pics/total_ep.encode.pdf', width=5, height=4, pointsize=12)
df=sort(sapply(ep_links,function(x) length(unlist(x))))
df1 = df
hist(df1,breaks=100,xlab='#EP links',ylab='#Tissues',main=NA);#axis(1,at=seq(0,2000,500),labels=c("0","500","1000","1500",">2000"))
abline(v=mean(df),col='red');abline(v=median(df),col='blue');legend('topright',legend=c("mean","median"),lty=c(1,1),lwd=c(2.5,2.5),col=c("red","blue"))
dev.off()


dist = sapply(1:length(ep_links),function(i){x=list.e.active.cell[[i]]; x[sapply(x,is.null)]=NULL; p_rep=rep(names(x),times=sapply(x,length)); abs(start(prom.bs.mid[p_rep])-start(enh.bs.mid[unlist(x)]))})
dist = unlist(dist)

## Fig. S2A - Histogram of the distances between E and P in all EP links across all cell types
pdf(file='pics/total_ep_dist.encode.pdf', width=5, height=4, pointsize=12)
df=dist
df[which(df>=1e5)]=1e5
hist(df,breaks=100,xlab='EP Distance [bp]',ylab='Frequency',main=NA,xaxt="n");axis(1,at=seq(2e4,1e5,2e4),labels=c("2e+04","4e+04","6e+04","8e+04",">1e+05"))
abline(v=mean(dist),col='red');abline(v=median(dist),col='blue');legend('topright',legend=c("mean","median"),lty=c(1,1),lwd=c(2.5,2.5),col=c("red","blue"))
dev.off()


