# Calculate for al samples belonging to the same cell type, s.cell, their distances from the other
# samples. Then it perform column means in case there are multiple samples from the same s.cell
calcDistFromSamples <- function(s.cell,y,...){
  dist.vals <- as.matrix(dist(t(y),diag=T)/sqrt(dim(y)[1]))
  cell.vals <- colMeans(dist.vals[rownames(dist.vals)==s.cell,,drop=F])
  cell.vals <- cell.vals[-which(names(cell.vals)==s.cell)]
  return(cell.vals)
}

# Calculate the specificity of the interactions predicted using the Shannon entropy (SE) 
# The lower the SE the more specific are the predicted EP links
shannon.entropy <- function(p)
{
  p <- p+abs(min(p))
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

# Create Dijk values - the multiplication between Ei and Pj in sample k, and returns Dij matrix
createDijMat <- function(ep.list,matE,matP,col2type,...){
  y <- lapply(names(ep.list),function(id){ x<-as.matrix(matE[ep.list[[id]],]); prod <- t(matrix(unlist(rep(matP[id,],dim(x)[1])),ncol=dim(x)[1]));
  res <- x*prod; return(res)})
  names(y) <- names(ep.list)
  x <- do.call(rbind,y)
  colnames(x) <- col2type[colnames(x),"col2type"]
  x = log2(x+1)
  return(x)
}

# Create Dijk values - the summation between Ei and Pj in sample k, and returns Dij matrix
createDijMatLog <- function(ep.list,matE,matP,col2type,...){
  y <- lapply(names(ep.list),function(id){ e.ids=ep.list[[id]]; if(any(is.na(e.ids))) e.ids=e.ids[!is.na(e.ids)]; x<-as.matrix(matE[e.ids,,drop=F]); mat=matP[id,,drop=F]; mat=mat[rep(seq_len(nrow(mat)), each=dim(x)[1]),]; res <- x + mat; return(res)})
  names(y) <- names(ep.list)
  x <- do.call(rbind,y)
  colnames(x) <- col2type[colnames(x),"cell_type"]
  return(x)
}

# Creates a set of TP EP links based on ChIA-PET and/or eQTL data
createTPepSet <- function(ep.list,g_to_e.chia.cell=NA,enh_name_hit_tmp=NA,...){
  true_ep <- ep.list
  eqtl_or_chia <- NA; only_eqtl=NA; only_chia=NA;
  if(!is.na(g_to_e.chia.cell)){
    true_ep <- lapply(names(ep.list), function(id){enh_ids = g_to_e.chia.cell[[id]]$eHit; if(is.na(enh_ids[1])) return(NULL); return(enh_ids[ep.list[[id]]]) })
    names(true_ep) <- names(ep.list)
  }
  if(!is.na(g_to_e.chia.cell) && !is.na(enh_name_hit_tmp)){
    eqtl_or_chia <- sapply(seq_len(length(true_ep)),function(i){res <- true_ep[[i]]; tmp=NULL; if(!is.null(enh_name_hit_tmp[[i]])) tmp=names(enh_name_hit_tmp[[i]])[enh_name_hit_tmp[[i]]==TRUE]; id=union(tmp,names(true_ep[[i]])[true_ep[[i]]==TRUE]); 
    res[id]=TRUE; return(res)})
    only_chia <- sapply(seq_len(length(true_ep)),function(i){res <- true_ep[[i]]; res[res==TRUE]=FALSE; tmp=NULL; if(!is.null(enh_name_hit_tmp[[i]])) tmp=names(enh_name_hit_tmp[[i]])[enh_name_hit_tmp[[i]]==TRUE]; id=setdiff(names(true_ep[[i]][true_ep[[i]]==TRUE]),tmp); 
    res[id]=TRUE; return(res)}) 
  }
  if(!is.na(enh_name_hit_tmp)){
    only_eqtl <- sapply(seq_len(length(true_ep)),function(i){res <- true_ep[[i]]; res[res==TRUE]=FALSE; tmp=NULL; if(!is.null(enh_name_hit_tmp[[i]])) tmp=names(enh_name_hit_tmp[[i]])[enh_name_hit_tmp[[i]]==TRUE]; id=setdiff(tmp,names(true_ep[[i]][true_ep[[i]]==TRUE])); 
    res[id]=TRUE; return(res)})
  }

  return(list(true_ep,eqtl_or_chia,only_eqtl,only_chia))
  
}

cellSpecificityScore <- function(s.cell, y, ...){
  dist = calcDistFromSamples(s.cell,y)
  dist = sapply(unique(names(dist)), function(cell){return(mean(dist[names(dist)==cell]))})
  tot_dist = sum(dist)
  x_c = rowMeans(y[,which(colnames(y)==s.cell),drop=F])
  if(dim(y)[1]!=1)
	x_i = sapply(unique(colnames(y)), function(cell) {ids=which(colnames(y)==cell); return(rowMeans(y[,ids,drop=F]))})
  else{
	x_i = sapply(unique(colnames(y)), function(cell) {ids=which(colnames(y)==cell); return(rowMeans(y[,ids,drop=F]))},USE.NAMES=FALSE)
	tf = rownames(y)
	x_i = matrix(x_i,nrow=1)
	rownames(x_i) = tf
	colnames(x_i) = unique(colnames(y))
  }
  x_i = x_i[,-which(colnames(x_i)==s.cell),drop=F]
  diff_x = (x_c-x_i)
  Sc=sum(t(diff_x)*dist)/tot_dist
  return(Sc)
}

range01 <- function(x){ (x-min(x))/(max(x)-min(x))}
