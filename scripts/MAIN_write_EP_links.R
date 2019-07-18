##################################################################################################################################
########################### Writing the predicted cell type-specific EP links as BED 12 format ###################################
##################################################################################################################################

require(biomaRt)

ensembl=useMart("ensembl")  # using ensembl database data
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)   # from ensembl using homosapien gene data

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

# Enhancer genomic positions
enh.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.enh.pos.rds',sep=''))
# Promoters positions
prom.bs <- readRDS(paste(data_directory,data_type,'/',data_type,'.prom.pos.rds',sep=''))

# Predicted cell type-specific EP links from CT-FOCS
ep_links = readRDS(file=paste0(data_directory,data_type,'/',data_type,'.EP_links.cell.rds'))

prom.bs.ep <- prom.bs[unlist(sapply(ep_links,function(x) rep(names(x),times=sapply(x,length))))]
prom.bs.ep$cell_type = unlist(sapply(names(ep_links),function(x) rep(x,length(unlist(ep_links[[x]])))))

enh.bs.ep <- enh.bs[unlist(sapply(ep_links,unlist))]

prom.bs.ep.mid = prom.bs.ep

if(data_type=='encode' | data_type=='roadmap'){
shift.size = abs(start(prom.bs.ep)-end(prom.bs.ep))/2
prom.bs.ep.mid = GenomicRanges::shift(prom.bs.ep,shift=shift.size)
prom.bs.ep.mid = promoters(prom.bs.ep.mid,upstream=0,downstream=1)
}

shift.size = abs(start(enh.bs.ep)-end(enh.bs.ep))/2
enh.bs.ep.mid = GenomicRanges::shift(enh.bs.ep,shift=shift.size)
enh.bs.ep.mid = promoters(enh.bs.ep.mid,upstream=0,downstream=1)

p_chr = as.character(seqnames(prom.bs.ep))
p_start = start(prom.bs.ep)
p_end = end(prom.bs.ep)
p_strand = rep('*',length(prom.bs.ep))
if(data_type=='fantom5')
	p_strand = strand(prom.bs.ep)

e_chr = as.character(seqnames(enh.bs.ep))
e_start = start(enh.bs.ep)
e_end = end(enh.bs.ep)
e_strand = rep('*',length(enh.bs.ep))

cell_type = prom.bs.ep$cell_type
prom_id = names(prom.bs.ep)
enh_id = names(enh.bs.ep)

ep_dist = start(prom.bs.ep.mid)-start(enh.bs.ep.mid)

if(data_type=='fantom5')
	ep_dist[which(strand(prom.bs.ep.mid)=='-')] = ep_dist[which(strand(prom.bs.ep.mid)=='-')]*(-1)

sym_id = sapply(prom.bs.ep$sym_id,function(x) paste0(x,collapse=';'))
ref_id = c()
if(data_type=='fantom5')
	ref_id = sapply(prom.bs.ep$ref_id,function(x) paste0(x,collapse=';'))
if(data_type=='encode' | data_type=='roadmap')
	ref_id = sapply(prom.bs.ep$ref_mrna_list,function(x) paste0(x,collapse=';'))

ent_id = sapply(prom.bs.ep$ent_id,function(x) paste0(x,collapse=';'))
ens_id = sapply(prom.bs.ep$ens_id,function(x) paste0(x,collapse=';'))


df = data.frame(p_chr=p_chr,p_start=p_start,p_end=p_end,p_strand=p_strand,e_chr=e_chr,e_start=e_start,e_end=e_end,e_strand=e_strand,cell_type=cell_type,prom_id=prom_id,enh_id=enh_id,
			ep_dist=ep_dist,sym_id=sym_id,ent_id=ent_id,ens_id=ens_id,ref_id=ref_id)

write.table(df, file=paste0(data_directory,data_type,'.cell_type_specific_ep_links.txt'),sep="\t",quote=FALSE,row.names=FALSE)



##################################### Do not Run ########################################

if(data_type=='fantom5'){
	gene_mrna = unique(unlist(prom.bs$ref_id))
	
	genes.with.id=getBM(attributes=c("entrezgene", "external_gene_name","ensembl_gene_id","refseq_mrna"),filters="refseq_mrna",
                    values=gene_mrna, mart= ensembl) # fuction to get  gene id's and gene name from data base
	res = mclapply(prom.bs$ref_id,function(x) {x=unique(unlist(x)); id=match(genes.with.id$refseq_mrna,x); id=which(!is.na(id)); if(length(id)==0) return(list(NULL,NULL)); 
			return(list(unique(genes.with.id[id,'entrezgene']),unique(genes.with.id[id,'ensembl_gene_id'])))},mc.cores=40)
	ent_id = sapply(res,function(x) x[1])
	ens_id = sapply(res,function(x) x[2])
	prom.bs$ent_id = ent_id
	prom.bs$ens_id = ens_id
	
}




