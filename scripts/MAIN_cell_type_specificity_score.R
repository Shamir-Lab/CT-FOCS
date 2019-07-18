# load objects and read files
data_directory = '../data/'
script_directory = '../scripts/'
tmp_directory = '../tmp/'
# data types - fantom5, encode
data_type <- 'encode'
EP_num_th <- 50 #Minimum number of cell type-specific EP links for identifying cell type-specific TFs 

## Specificity source functions
source(paste0(script_directory,'FUNCTIONS_specificity_score.R'))

#############################################################################################
## EP signal - ENCODE
#############################################################################################


exp.data <- readRDS(paste0(data_directory,'Sheffield_GR_2013.exp.data.rds'))
list.e.active.cell = readRDS(file=paste0(data_directory,data_type,'.list.e.active.cell.rds'))

cells_to_check <- names(list.e.active.cell)[sapply(list.e.active.cell,function(x) length(unlist(x))>=EP_num_th)]
Sc_ranks = rep(length(list.e.active.cell),length(cells_to_check))
names(Sc_ranks) = cells_to_check
ep_dist_per_cell = rep(list(NA),length(cells_to_check))
names(ep_dist_per_cell) = cells_to_check

for( i in 1:length(cells_to_check)){
	print(i)
	s <- cells_to_check[i]

	e.active.cell <- list.e.active.cell[[s]]
	p.active.cell <- names(e.active.cell)
	length(p.active.cell)
	
	tmp <- e.active.cell
	tmp[sapply(tmp,is.null)] <- NULL
	length(unlist(tmp)) # num of EP
	length(unique(unlist(tmp))) # num of E
	length(tmp) # num of P
	
	ep.list <- tmp
	
	
	prom.ep <- start(prom.bs.mid[rep(names(tmp),times=sapply(tmp,length))])
	enh.ep <- start(enh.bs.mid[unlist(tmp)])
	ep_dist_per_cell[[s]] = abs(enh.ep-prom.ep)

}

pdf(file='pics/ep_dist.encode.rle.by.pdf', width=5, height=4, pointsize=12)
hist(unlist(ep_dist_per_cell),breaks=100, xlab="EP distance [bp]", main="ENCODE EP distances")
dev.off()

#############################################################################################
## GE - ENCODE
#############################################################################################

colnames(exp.data) = gsub("\\.","-",colnames(exp.data))

Sc_exp_ranks = rep(length(list.e.active.cell),length(cells_to_check))
names(Sc_exp_ranks) = cells_to_check

for( i in 1:length(cells_to_check)){
	print(i)
	s <- cells_to_check[i]

	e.active.cell <- list.e.active.cell[[s]]
	p.active.cell <- names(e.active.cell)
	length(p.active.cell)
	
	tmp <- e.active.cell
	tmp[sapply(tmp,is.null)] <- NULL
	length(unlist(tmp)) # num of EP
	length(unique(unlist(tmp))) # num of E
	length(tmp) # num of P
	
	ep.list <- tmp

	
	length(gene.sym <- unique(unlist(prom.bs[names(tmp)]$sym_id)))
	exp.data.tmp <- exp.data[which(!is.na(match(rownames(exp.data),gene.sym))),]

	y = as.matrix(exp.data.tmp)
	
	Sc = mclapply(unique(colnames(y)), function(cell) cellSpecificityScore(cell,y), mc.cores=10)
	names(Sc) = unique(colnames(y))
	Sc = unlist(Sc)
	Sc = sort(Sc,decreasing=T)
	Sc.col = rep("black",length(Sc))
	names(Sc.col) <- names(Sc)
	
	Sc_exp_ranks[s] = which(names(Sc)==s)
	print(paste0(s, " rank: ", Sc_exp_ranks[s]))
}

pdf("pics/encode.cellType.spec.score.ranks.pdf", width=10, height=5)
par(mfrow=c(1,2))
barplot(table(Sc_ranks),xlab="EP signal Sc rank (out of 106 cell types)", ylab="#cell-types")
barplot(table(Sc_exp_ranks),xlab="GE Sc rank (out of 112 cell types)", ylab="#cell-types")
dev.off()


saveRDS(Sc_ranks,file="RDS/encode.ep.sc.scores.rds")
saveRDS(Sc_exp_ranks,file="RDS/encode.ge.sc.scores.rds")
saveRDS(ep_dist_per_cell,file="RDS/ep_dist_per_cell.encode.rds")


#############################################################################################
## FANTOM5 CT-FOCS - EP signal
#############################################################################################
data_type <- 'fantom5'
list.e.active.cell = readRDS(file=paste0(data_directory,data_type,'.list.e.active.cell.rds'))

cells_to_check <- names(list.e.active.cell)[sapply(list.e.active.cell,function(x) length(unlist(x))>=EP_num_th)]
Sc_ranks = rep(length(list.e.active.cell),length(cells_to_check))
names(Sc_ranks) = cells_to_check
ep_dist_per_cell = rep(list(NA),length(cells_to_check))
names(ep_dist_per_cell) = cells_to_check

for( i in 1:length(cells_to_check)){
	if(i%%10==0) print(i)
	s <- cells_to_check[i]

	e.active.cell <- list.e.active.cell[[s]]
	p.active.cell <- names(e.active.cell)
	length(p.active.cell)
	
	tmp <- e.active.cell
	tmp[sapply(tmp,is.null)] <- NULL
	tmp <- sapply(tmp,function(x){id=which(is.na(x)); if(length(id)!=0) return(x[!id]) else return(x)})
	length(unlist(tmp)) # num of EP
	length(unique(unlist(tmp))) # num of E
	length(tmp) # num of P
	
	ep.list.tmp <- tmp

	prom.ep <- start(prom.bs.mid[rep(names(ep.list.tmp),times=sapply(tmp,length))])
	enh.ep <- start(enh.bs.mid[unlist(ep.list.tmp)])
	ep_dist_per_cell[[s]] = abs(enh.ep-prom.ep)
	
	
	genes <- intersect(names(ep.list.tmp),names(ep.list))
	if(length(genes)==0) next
	res1 = sapply(genes,function(x) length(setdiff(ep.list[[x]],ep.list.tmp[[x]])))
	res2 = sapply(genes,function(x) length(setdiff(ep.list.tmp[[x]],ep.list[[x]])))
	res3 = res1>0 | res2>0
	num_diff_E[[s]] = list(res1,res2,res3)
}

#############################################################################################
## GE CT-FOCS
#############################################################################################

num_ep_th = 50
colnames(exp.data) = tolower(gsub("\\.","-",colnames(exp.data)))
colnames(exp.data)[85] = "mcf7"

Sc_exp_ranks = rep(length(list.e.active.cell),length(cells_to_check))
names(Sc_exp_ranks) = cells_to_check

for( i in 1:length(cells_to_check)){
	print(i)
	s <- cells_to_check[i]

	e.active.cell <- list.e.active.cell[[s]]
	p.active.cell <- names(e.active.cell)
	length(p.active.cell)
	
	tmp <- e.active.cell
	tmp[sapply(tmp,is.null)] <- NULL
	length(unlist(tmp)) # num of EP
	length(unique(unlist(tmp))) # num of E
	length(tmp) # num of P
	
	ep.list <- tmp

	
	length(gene.sym <- unique(unlist(prom.bs[names(tmp)]$sym_id)))
	exp.data.tmp <- exp.data[which(!is.na(match(rownames(exp.data),gene.sym))),]

	y = as.matrix(exp.data.tmp)
	
	Sc = mclapply(unique(colnames(y)), function(cell) cellSpecificityScore(cell,y), mc.cores=10)
	names(Sc) = unique(colnames(y))
	Sc = unlist(Sc)
	Sc = sort(Sc,decreasing=T)
	Sc.col = rep("black",length(Sc))
	names(Sc.col) <- names(Sc)
	
	Sc_exp_ranks[s] = which(names(Sc)==s)
	print(paste0(s, " rank: ", Sc_exp_ranks[s]))
}


saveRDS(Sc_ranks,file="RDS/fantom.ep.sc.scores.ct_focs.rds")
saveRDS(Sc_exp_ranks,file="RDS/fantom.ge.sc.scores.ct_focs.rds")

#############################################################################################
## FANTOM5 JEME - GE
#############################################################################################

colnames(exp.data) = tolower(gsub("\\.","-",colnames(exp.data)))
colnames(exp.data)[85] = "mcf7"

score <- 0.3
tb = table(prom.gr.jeme$tissue)
cells_to_check <- names(tb)[sapply(tb,function(x) x>=EP_num_th)]
cells_to_check = intersect(cells_to_check,colnames(exp.data))
Sc_exp_ranks_jeme = rep(length(list.e.active.cell),length(cells_to_check))
names(Sc_exp_ranks_jeme) = cells_to_check

for( i in 1:length(cells_to_check)){
	print(i)
	s = cells_to_check[i]
	enh.gr.jeme.cell <- enh.gr.jeme.mid[unlist(mclapply(tolower(enh.gr.jeme$tissue),function(x) (tolower(s) %in% x),mc.cores=10))]
	prom.gr.jeme.cell <- prom.gr.jeme.mid[unlist(mclapply(tolower(prom.gr.jeme$tissue),function(x) (tolower(s) %in% x),mc.cores=10))]

	match.id <- which(prom.gr.jeme.cell$score>=score)
	enh.gr.jeme.cell <- enh.gr.jeme.cell[match.id]
	prom.gr.jeme.cell <- prom.gr.jeme.cell[match.id]

	ep.int <- paste0(names(enh.gr.jeme.cell),'##',names(prom.gr.jeme.cell))
	ep.int <- unique(ep.int)
	e.name <- unlist(strsplit(ep.int,"##"))[seq(1,2*length(ep.int),2)]
	p.name <- unlist(strsplit(ep.int,"##"))[seq(2,2*length(ep.int),2)]
	enh.gr.jeme.cell <- enh.gr.jeme.cell[e.name]
	prom.gr.jeme.cell <- prom.gr.jeme.cell[p.name]
	prom_enh.df.jeme <- data.frame(p.names=names(prom.gr.jeme.cell),e.names=names(enh.gr.jeme.cell))
	prom_enh.list.jeme <- dlply(prom_enh.df.jeme, .(p.names), function(dt) return(as.character(dt$e.names)))
	
	ep.list <- prom_enh.list.jeme

	p.names <- unique(names(ep.list))
	gene.sym.jeme <- unique(unlist(sapply(p.names, function(x){ y<-unlist(strsplit(x,"%")); y<-unlist(strsplit(y,"\\$")); return(y[seq(2,length(y),5)])})))

	length(gene.sym <- unique(unlist(prom.bs[names(tmp)]$sym_id)))
	exp.data.tmp <- exp.data[which(!is.na(match(rownames(exp.data),gene.sym.jeme))),]

	y = as.matrix(exp.data.tmp)
	
	Sc = mclapply(unique(colnames(y)), function(cell) cellSpecificityScore(cell,y), mc.cores=10)
	names(Sc) = unique(colnames(y))
	Sc = unlist(Sc)
	Sc = sort(Sc,decreasing=T)
	Sc.col = rep("black",length(Sc))
	names(Sc.col) <- names(Sc)
	
	Sc_exp_ranks_jeme[s] = which(names(Sc)==s)
	print(paste0(s, " rank: ", Sc_exp_ranks_jeme[s]))
}

#############################################################################################
## FANTOM5 JEME - EP signal
#############################################################################################

score <- 0.3
tb = table(prom.gr.jeme$tissue)
Sc_ranks_jeme = rep(length(list.e.active.cell),length(cells_to_check))
names(Sc_ranks_jeme) = cells_to_check

for( i in 1:length(cells_to_check)){
	if(i%%10==0) print(i)
	s <- cells_to_check[i]

	enh.gr.jeme.cell <- enh.gr.jeme.mid[unlist(mclapply(tolower(enh.gr.jeme$tissue),function(x) (tolower(s) %in% x),mc.cores=10))]
	prom.gr.jeme.cell <- prom.gr.jeme.mid[unlist(mclapply(tolower(prom.gr.jeme$tissue),function(x) (tolower(s) %in% x),mc.cores=10))]

	match.id <- which(prom.gr.jeme.cell$score>=score)
	enh.gr.jeme.cell <- enh.gr.jeme.cell[match.id]
	prom.gr.jeme.cell <- prom.gr.jeme.cell[match.id]

	ep.int <- paste0(names(enh.gr.jeme.cell),'##',names(prom.gr.jeme.cell))
	ep.int <- unique(ep.int)
	e.name <- unlist(strsplit(ep.int,"##"))[seq(1,2*length(ep.int),2)]
	p.name <- unlist(strsplit(ep.int,"##"))[seq(2,2*length(ep.int),2)]
	enh.gr.jeme.cell <- enh.gr.jeme.cell[e.name]
	prom.gr.jeme.cell <- prom.gr.jeme.cell[p.name]
	prom_enh.df.jeme <- data.frame(p.names=names(prom.gr.jeme.cell),e.names=names(enh.gr.jeme.cell))
	prom_enh.list.jeme <- dlply(prom_enh.df.jeme, .(p.names), function(dt) return(as.character(dt$e.names)))
	
	ep.list <- prom_enh.list.jeme

	x <- createDijMatLog(ep.list,Me,Mg,col2type)
	y <- x


	Sc = mclapply(unique(colnames(y)), function(cell) cellSpecificityScore(cell,y), mc.cores=10)
	names(Sc) = unique(colnames(y))
	Sc = unlist(Sc)
	Sc = sort(Sc,decreasing=T)
	Sc.col = rep("black",length(Sc))
	names(Sc.col) <- names(Sc)
	
	Sc_ranks_jeme[s] = which(names(Sc)==s)
	print(paste0(s, " rank: ", Sc_ranks_jeme[s]))
}

saveRDS(Sc_ranks_jeme,file="RDS/fantom.ep.sc.scores.jeme.rds")
saveRDS(Sc_exp_ranks_jeme,file="RDS/fantom.ge.sc.scores.jeme.rds")

Sc_ranks_jeme=readRDS(file="RDS/fantom.ep.sc.scores.jeme.rds")
Sc_ranks=readRDS(file="RDS/fantom.ep.sc.scores.ct_focs.rds")

Sc_exp_ranks=readRDS(file="RDS/fantom.ge.sc.scores.ct_focs.rds")
Sc_exp_ranks_jeme=readRDS(file="RDS/fantom.ge.sc.scores.jeme.rds")

cells = intersect(names(Sc_ranks_jeme),names(Sc_ranks))

df_ep = data.frame("CT_FOCS"=Sc_ranks[cells],'JEME'=Sc_ranks_jeme[cells])
head(df_ep[order(df_ep$CT_FOCS),],20)

df_exp = data.frame("CT_FOCS"=Sc_exp_ranks[cells],'JEME'=Sc_exp_ranks_jeme[cells])
head(df_exp[order(df_exp$CT_FOCS),],20)


df = data.frame('Sc_rank'=c(Sc_ranks[cells],Sc_ranks_jeme[cells]),'Method'=c(rep('CT-FOCS',length(cells)),rep('JEME',length(cells))))
cdf = data.frame(Method=c('CT-FOCS','JEME'), Sc_mean=c(mean(Sc_ranks[cells]),mean(Sc_ranks_jeme[cells])))


pdf("pics/fantom.cellType.spec.score.ranks.ct_focs_jeme.density.pdf", width=10, height=5)
p1 = ggplot(data=df,aes(x=Sc_rank,fill=Method, colour=Method))+ geom_density(alpha=0.3)+ geom_vline(data=cdf, aes(xintercept=Sc_mean,  colour=Method),
               linetype="dashed", size=1)+ xlab("Cell Type Specificity Score Rank (out of 328 cell types)") + scale_fill_manual(values=c("#1a85ff","#d41159")) + scale_colour_manual(values=c("#1a85ff","#d41159"))# + theme(legend.pos = "none")

#p2 = ggplot(data=df,aes(y=Sc_rank,x=Method, fill=Method, colour=Method))+ geom_boxplot(alpha=0.5, outlier.shape=16, outlier.size=1, notch=TRUE) + ylab("Cell Type Specificity Score Rank (out of 328 cell types)")
#p2 = p2 + geom_jitter(shape=16, position=position_jitter(0.2)) + scale_fill_manual(values=c("#1a85ff","#d41159")) + scale_colour_manual(values=c("#1a85ff","#d41159"))

p1 + theme(axis.text.x = element_text( color="black", 
                           size=14),
          axis.text.y = element_text( color="black", 
                           size=14))
#multiplot(p1, p2, cols=2)
dev.off()

p.v=wilcox.test(Sc_ranks[cells],Sc_ranks_jeme[cells],paired=TRUE,alternative="less")$p.value

pdf("pics/fantom.cellType.spec.boxplot.score.ranks.pdf")
boxplot(list('CT-FOCS'=Sc_ranks[cells],JEME=Sc_ranks_jeme[cells]),ylab="Cell Type Specificity Score Rank (out of 328 cell types)")
text(coords2,labels="P<9.6E-39")
dev.off()

coords2 = locator(1)



rel_specif_cell_jeme = readRDS(file=paste0('RDS/rel_specif_cell_jeme.',s.cell,'.rds'))
rel_specif_cell = readRDS(file=paste0('RDS/rel_specif_cell.',s.cell,'.rds'))


cells = intersect(names(rel_specif_cell),names(rel_specif_cell_jeme))

df = data.frame('Rel_Exp'=log2(c(rel_specif_cell[cells],rel_specif_cell_jeme[cells])),'Method'=c(rep('CT-FOCS',length(cells)),rep('JEME',length(cells))))
cdf = data.frame(Method=c('CT-FOCS','JEME'), Rel_exp_mean=c(mean(log2(rel_specif_cell[cells])),mean(log2(rel_specif_cell_jeme[cells]))))

p.v=wilcox.test(rel_specif_cell[cells],rel_specif_cell_jeme[cells],paired=TRUE,alternative='greater')$p.value

pdf(paste0("pics/",s.cell,".relative_precision_fantom.ctfocs.jeme.pdf"), width=10, height=5)

p2 = ggplot(data=df,aes(y=Rel_Exp,x=Method, fill=Method, colour=Method))+ geom_boxplot(alpha=0.5, outlier.shape=16, outlier.size=1, notch=TRUE) + ylab("GM12878 Enrichment [log2[fold change]]")
p2 = p2 + geom_jitter(shape=16, position=position_jitter(0.2)) + scale_fill_manual(values=c("#1a85ff","#d41159")) + scale_colour_manual(values=c("#1a85ff","#d41159"))

dev.off()

