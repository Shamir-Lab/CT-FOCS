##################################################################################################################################
################################# External validation using ChIA-PET connected loop sets #########################################
##################################################################################################################################

########################################## Note: Only works via R unix platform ##################################################


## This is a line by line running tutorial. In order to to run a subset of lines please select them using shift and run using ctrl+R

## Please follow the instruction given in http://acgt.cs.tau.ac.il/ct-focs/tutorial.html

# load libraries - make sure you pre-installed these packages
cran_libs = c('S4Vectors','MASS','parallel','igraph','MatchIt','RColorBrewer','ggplot2')
bioconductor_libs = c('GenomicRanges')
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

source(paste0(script_directory,'FUNCTIONS_ChIA_PET_external_validation.R'))

# Enhancer genomic positions
enh.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.enh.pos.rds',sep=''))
# Promoters positions
prom.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.prom.pos.rds',sep=''))

# Promoter to closest enhnacers in fixed window:
g_to_e = readRDS(paste(data_directory,data_type,'/',data_type,'.cand.enh.rds',sep=''))

# Predicted cell type-specific EP links from CT-FOCS
ep_links = readRDS(file=paste0(data_directory,data_type,'/',data_type,'.EP_links.cell.rds'))

# load processed ChIA-PET loops
chia.1.bs <- readRDS(file=paste0(data_directory,"chia_pet_specific.1.noOverlaps.anchor.5kb.rds"))
chia.2.bs <- readRDS(file=paste0(data_directory,"chia_pet_specific.2.noOverlaps.anchor.5kb.rds"))
unique(chia.1.bs$cellType)

# load GM12878 topological associated domains (TADs)
df <- read.table(paste0(data_directory,'osborne_tads_gm12878.txt'),header=TRUE)
hic.bs <- GRanges(seqnames=df$Chr, ranges=IRanges(start=df$start,end=df$end),strand=rep('*',dim(df)[1]))

# we work on GM12878 ChIA-PET loops
s.cell <- 'GM12878'
chia.1.bs.cell <- chia.1.bs[chia.1.bs$cellType==s.cell]
chia.2.bs.cell <- chia.2.bs[chia.2.bs$cellType==s.cell]

shift.size=abs(start(chia.1.bs.cell)-end(chia.1.bs.cell))/2
chia.1.bs.cell.mid = shift(chia.1.bs.cell,shift=shift.size)
start.pos = start(chia.1.bs.cell.mid)
chia.bs.cell.mid = promoters(chia.1.bs.cell.mid,upstream=0,downstream=1)
chia.1.bs.cell.mid = promoters(chia.1.bs.cell.mid,upstream=win.size,downstream=win.size)


shift.size=abs(start(chia.2.bs.cell)-end(chia.2.bs.cell))/2
chia.2.bs.cell.mid = shift(chia.2.bs.cell,shift=shift.size)
end.pos = start(chia.2.bs.cell.mid)
chia.bs.cell.mid = c(chia.bs.cell.mid,promoters(chia.2.bs.cell.mid,upstream=0,downstream=1))
chia.2.bs.cell.mid = promoters(chia.2.bs.cell.mid,upstream=win.size,downstream=win.size)


orig.dis = abs(start.pos-end.pos)


names(chia.1.bs.cell.mid)=paste0(seqnames(chia.1.bs.cell),':',start(chia.1.bs.cell),'-',end(chia.1.bs.cell))
names(chia.2.bs.cell.mid)=paste0(seqnames(chia.2.bs.cell),':',start(chia.2.bs.cell),'-',end(chia.2.bs.cell))
names(chia.bs.cell.mid) = c(names(chia.1.bs.cell.mid),names(chia.2.bs.cell.mid))
chia.bs.cell.mid = c(chia.1.bs.cell.mid,chia.2.bs.cell.mid)

shift.size=abs(start(enh.bs)-end(enh.bs))/2
enh.bs.mid = shift(enh.bs,shift=shift.size)
enh.bs.mid = promoters(enh.bs.mid,upstream=win.size,downstream=win.size)

shift.size=abs(start(prom.bs)-end(prom.bs))/2
prom.bs.mid = shift(prom.bs,shift=shift.size)
prom.bs.mid = promoters(prom.bs.mid,upstream=win.size,downstream=win.size)

#remove chia loops crossing TADs
hit.1 = findOverlaps(chia.1.bs.cell.mid,hic.bs,ignore.strand=TRUE)
hit.2 = findOverlaps(chia.2.bs.cell.mid,hic.bs,ignore.strand=TRUE)

not.in.hic.bs <- 1:length(chia.1.bs.cell.mid)
not.in.hic.bs <- setdiff(not.in.hic.bs,union(unique(queryHits(hit.1)),unique(queryHits(hit.2))))

s.hit.1 <- unique(subjectHits(hit.1))
s.hit.2 <- unique(subjectHits(hit.2))
s.hit <- union(s.hit.1,s.hit.2)

intra.loops = lapply(1:length(s.hit),function(i){if(i%%100==0) print(i); s=s.hit[i]; return(intersect(queryHits(hit.1)[subjectHits(hit.1)==s],queryHits(hit.2)[subjectHits(hit.2)==s]))})
intra.loops = unique(do.call(c,intra.loops))
length(intra.loops)/length(chia.1.bs.cell.mid)
intra.loops <- unique(c(intra.loops,not.in.hic.bs))

chia.1.bs.cell.mid = chia.1.bs.cell.mid[intra.loops]
chia.2.bs.cell.mid = chia.2.bs.cell.mid[intra.loops]
chia.bs.cell.mid = c(chia.1.bs.cell.mid,chia.2.bs.cell.mid)


## create bipartite graph of the chiapet interactions
types = rep(c(0,1),times=length(chia.1.bs.cell.mid))
edges = seq(1,2*length(chia.1.bs.cell.mid),by=1)
chia.gp = make_bipartite_graph(types, edges, directed = FALSE)
vertex.name.0 = names(chia.1.bs.cell.mid)#paste0(seqnames(chia.1.bs.cell),':',start(chia.1.bs.cell),'-',end(chia.1.bs.cell))
vertex.name.1 = names(chia.2.bs.cell.mid)#paste0(seqnames(chia.2.bs.cell),':',start(chia.2.bs.cell),'-',end(chia.2.bs.cell))
vertex.name = rep(NA,2*length(chia.1.bs.cell.mid))
vertex.name[types==0] = vertex.name.0
vertex.name[types==1] = vertex.name.1

chia.gp = set.vertex.attribute(chia.gp, "name", value=vertex.name)

hit.1.1 <- findOverlaps(chia.1.bs.cell.mid,chia.1.bs.cell.mid,ignore.strand=TRUE,minoverlap=win.size/2)
hit.1.2 <- findOverlaps(chia.1.bs.cell.mid,chia.2.bs.cell.mid,ignore.strand=TRUE,minoverlap=win.size/2)
hit.2.2 <- findOverlaps(chia.2.bs.cell.mid,chia.2.bs.cell.mid,ignore.strand=TRUE,minoverlap=win.size/2)

#add edges to the graph - overlapping anchros
# if we find an overlap b/w different anchors than we add an edge b/w one overlapping respective anchors from chia.2.bs.cell.mid 
df = as.data.frame(hit.1.1)
keep <- apply(df, 1, function(x) length(unique(x[!is.na(x)])) != 1)
df = df[keep,]
df1 = df
df2 = df
df1[,1] = names(chia.1.bs.cell.mid)[df[,1]]
df1[,2] = names(chia.2.bs.cell.mid)[df[,2]]
df2[,1] = names(chia.2.bs.cell.mid)[df[,1]]
df2[,2] = names(chia.1.bs.cell.mid)[df[,2]]
df <- rbind(df1,df2)
df.flat = unlist(apply(df,1,function(x) list(c(x[1],x[2]))))
chia.gp <- add_edges(chia.gp, edges=df.flat)

df = as.data.frame(hit.1.2)
keep <- apply(df, 1, function(x) length(unique(x[!is.na(x)])) != 1)
df = df[keep,]
df1 = df
df2 = df
df1[,1] = names(chia.1.bs.cell.mid)[df[,1]]
df1[,2] = names(chia.1.bs.cell.mid)[df[,2]]
df2[,1] = names(chia.2.bs.cell.mid)[df[,1]]
df2[,2] = names(chia.2.bs.cell.mid)[df[,2]]
df <- rbind(df1,df2)
df.flat = unlist(apply(df,1,function(x) list(c(x[1],x[2]))))
chia.gp <- add_edges(chia.gp, edges=df.flat)

df = as.data.frame(hit.2.2)
keep <- apply(df, 1, function(x) length(unique(x[!is.na(x)])) != 1)
df = df[keep,]
df1 = df
df2 = df
df1[,1] = names(chia.2.bs.cell.mid)[df[,1]]
df1[,2] = names(chia.1.bs.cell.mid)[df[,2]]
df2[,1] = names(chia.1.bs.cell.mid)[df[,1]]
df2[,2] = names(chia.2.bs.cell.mid)[df[,2]]
df <- rbind(df1,df2)
df.flat = unlist(apply(df,1,function(x) list(c(x[1],x[2]))))
chia.gp <- add_edges(chia.gp, edges=df.flat)

chrs = as.character(unique(seqnames(chia.1.bs.cell.mid)))

chia.gp = simplify(chia.gp)

comp = components(chia.gp)

saveRDS(comp$membership,file=paste0(data_directory,data_type,'.',s.cell,'.comp_list.',win.size*2,'.rds'))

## Go over the candidate list of enhancers per promoter and mark true positive EP links

g_list = names(g_to_e)
g_list = strsplit(g_list,"[_:]")
g_list = do.call(c,lapply(g_list,function(x) unlist(x)[2]))

comp_list = readRDS(file=paste0(data_directory,data_type,'.',s.cell,'.comp_list.',win.size*2,'.rds'))

for(chr in chrs){
	if(file.exists(paste0(tmp_directory,data_type,'.g_to_e.chia.cell.comp.',s.cell,'.',chr,'.rds'))) next;
	tmp_g_to_e <- g_to_e[which(g_list==chr)]
	print(paste0("Now working on: ",chr," ### total_promoters: ", length(tmp_g_to_e)))
	g_to_e.chia.cell.focs.path <- mclapply(seq_len(length(tmp_g_to_e)),function(i) make_g_to_e_ov_chia_pet_by_comp(i,s.cell,p_chr=names(tmp_g_to_e[i]),e_chrs=tmp_g_to_e[[i]],
							prom.bs.mid,enh.bs.mid,chia.bs.cell.mid,comp_list), mc.cores=30)
	names(g_to_e.chia.cell.focs.path) <- names(tmp_g_to_e)
	saveRDS(g_to_e.chia.cell.focs.path, file=paste0(tmp_directory,data_type,'.g_to_e.chia.cell.comp.',s.cell,'.',chr,'.rds'))
}

g_to_e.chia.cell.focs.path <- list()

for(chr in chrs){
	g_to_e.chia.cell.focs.path <- c(g_to_e.chia.cell.focs.path,	readRDS(file=paste0(tmp_directory,data_type,'.g_to_e.chia.cell.comp.',s.cell,'.',chr,'.rds')))
}

sum(sapply(g_to_e.chia.cell.focs.path,function(x) sum(x$eHit)))

saveRDS(g_to_e.chia.cell.focs.path, file=paste0(data_directory,data_type,'.chia_pet.connected.loop.sets.',s.cell,'.rds'))

g_to_e.chia.cell.focs.path = readRDS(file=paste0(data_directory,data_type,'.chia_pet.connected.loop.sets.',s.cell,'.rds'))

## find true positive EP links supported by ChIA-PET connected loop sets (CLSs)
s.cell <- 'GM12878'

e.active.cell <- ep_links[[s.cell]]
if(data_type=='fantom5')
	e.active.cell <- ep_links[[tolower(s.cell)]]

p.active.cell <- names(e.active.cell)
length(p.active.cell)

tmp <- e.active.cell
tmp[sapply(tmp,is.null)] <- NULL

length(unlist(tmp)) # num of EP links in s.cell
length(unique(unlist(tmp))) # num of linked enhancers in s.cell
length(tmp) # num of linked promoters in s.cell

ep.list = tmp

true_ep <- lapply(names(ep.list), function(id){ enh_ids = g_to_e.chia.cell.focs.path[[id]]$eHit; if(is.null(enh_ids)) return(NULL); if(is.na(enh_ids[1])) return(NULL); return(enh_ids[ep.list[[id]]]) })
names(true_ep) <- names(ep.list)
which.max(sapply(true_ep,sum))
Nt <- sum(unlist(true_ep)) #number of true positives - EP link supported by a CLS

#random permutations on chia_pet paths to see if we get significant EP pairs supported by chia_pet paths
set.seed(1)
perms = 1000
Nr <- rep(0,perms)

g_ep_list = names(ep.list)
g_ep_list = strsplit(g_ep_list,"[_:]")
g_ep_list = do.call(c,lapply(g_ep_list,function(x) unlist(x)[2]))

tmp_ep_list <- rep(list(NULL),length(chrs))
target_ep_links <- rep(list(NULL),length(chrs))
Tr <- rep(list(NULL),length(chrs))
names(tmp_ep_list) <- chrs
names(target_ep_links) <- chrs
names(Tr) <- chrs

for(chr in chrs){
	print(paste0("Now working on: ",chr))
	tmp_ep_list[[chr]] <- ep.list[which(g_ep_list==chr)]
	if(length(tmp_ep_list[[chr]])==0) next
	target_ep_links[[chr]] <- sapply(names(tmp_ep_list[[chr]]), function(p_chr) paste0(p_chr,'##',tmp_ep_list[[chr]][[p_chr]]))
	Tr[[chr]] <- createMatchItTargetVector(ep_dist.list[[chr]],unlist(target_ep_links[[chr]]))
}

summary(sapply(Tr,length))

set.seed(1)
for(i in 1:perms){
	print(paste0("Now working on permutation: ",i))
	for(chr in chrs){
		#print(paste0("Now working on: ",chr))
		if(length(tmp_ep_list[[chr]])==0) next
		tmp_ep_list.crt <- createMatchItCrtSet(ep_dist.list[[chr]], Tr[[chr]])
		tmp_ep_list.crt <- do.call(rbind,strsplit(tmp_ep_list.crt,"##"))
		true_ep <- apply(tmp_ep_list.crt,1,function(row) {enh_ids = g_to_e.chia.cell.focs.path[[row[1]]]$eHit; if(is.na(enh_ids[1])) return(NULL); return(enh_ids[row[2]])})
		Nr[i] <- Nr[i]+sum(true_ep)
	}
	print(Nr[i])
}

summary(Nr)

saveRDS(Nr, file=paste0(tmp_directory,data_type,".Nr.rand.perm.",s.cell,".rds"))

## Fig. S3B - histogram of random permutation results 
pdf(file=paste0('pics/',data_type,'.',s.cell,'.Nr.rand.perm.loops.pdf'))
hist(Nr,breaks=50,xlab="#EP supported by CLSs",main="Frequency of #EP supported by CLSs",xlim=c(0,Nt+50));abline(v=Nt,col='red');legend('top',legend="true value",col="red",lwd=1.5)
dev.off()


