### Functions - External validation with ChIA-PET connected loop sets ###

require(MatchIt)

###############################################################################################################################################
# The function checks for each candidate EP link if its E and P regions overlap with two different anchors of the same connected loop set (CLS)
###############################################################################################################################################
make_g_to_e_ov_chia_pet_by_comp <- function(i,s.cell,p_chr,e_chrs,prom.bs.mid,enh.bs.mid,chia.bs.cell.mid,comp_list,win.size=0.5e3, ...){
	if(i%%100==0) print(paste0('################### i: ',i))
	res <- rep(FALSE,length(e_chrs))
	names(res) <- e_chrs
	
	hit.p = findOverlaps(prom.bs.mid[p_chr],chia.bs.cell.mid,ignore.strand=TRUE,minoverlap=win.size/2)
	hit.e = findOverlaps(enh.bs.mid[e_chrs],chia.bs.cell.mid,ignore.strand=TRUE,minoverlap=win.size/2)
	
	if(dim(as.data.frame(hit.p))[1]==0) return(list(eHit=res))
	if(dim(as.data.frame(hit.e))[1]==0) return(list(eHit=res))
	
	s.hit.p = names(chia.bs.cell.mid)[unique(subjectHits(hit.p))]
	
	q.hit.e = unique(queryHits(hit.e))
	
	df.hit.e = data.frame(e=e_chrs[queryHits(hit.e)],s.hit.e=names(chia.bs.cell.mid)[subjectHits(hit.e)])
	flags = ddply(df.hit.e, .(e), .fun=function(dt){
		s.hit.e = dt$s.hit.e
		tmp_p_hit = s.hit.p
		tmp_e_hit = s.hit.e
		s.hit.e = setdiff(tmp_e_hit,tmp_p_hit)
		tmp_p_hit = setdiff(tmp_p_hit,tmp_e_hit)
	
		if(length(s.hit.e)==0 | length(tmp_p_hit)==0) return(FALSE)
		
		is.in.comp = length(intersect(comp_list[s.hit.e],comp_list[tmp_p_hit]))!=0
		return(c(flag=is.in.comp))
	})
	res[as.character(flags[,1])] <- flags[,2]
	return(list(eHit=res))
}

###############################################################################################################################################
# The function performs a single permutation using matchIt to match EP links by distance to random EP links with the same distance
# Then the random EP links are overlapped with ChIA-PET interactions 
# The function retuns a logical vector of EP links with TRUE as supported by 3D interactions, and FALSE otherwise
###############################################################################################################################################
make_single_random_perm_ov_chia_pet_by_win <- function(i,s.cell,p_chr,e_chrs,prom.bs.mid,enh.bs.mid,chia.1.bs.cell.mid,chia.2.bs.cell.mid,win.size=0.5e3,dists=NA,...){
	chr = unlist(strsplit(p_chr,"[_:]"))[2]
	enh_names <- names(enh.bs.mid)[which(seqnames(enh.bs.mid)==chr)]
	e_chrs_rand <- sample(x=enh_names,size=length(e_chrs))
	res <- make_g_to_e_ov_chia_pet_by_win(i,s.cell,p_chr=p_chr,e_chrs=e_chrs_rand,
							prom.bs.mid,enh.bs.mid,chia.1.bs.cell.mid,chia.2.bs.cell.mid,dists=dists)
}

###############################################################################################################################################
# The function performs multiple permutations (according to perms variable) using matchIt to match EP links 
# by distance to random EP links with the same distance
# Then the random EP links are overlapped with ChIA-PET interactions 
# The function retuns a logical vector of EP links with TRUE as supported by 3D interactions, and FALSE otherwise
###############################################################################################################################################
make_mult_random_perm_ov_chia_pet_by_win <- function(i,s.cell,p_chr,e_chrs,prom.bs.mid,enh.bs.mid,chia.1.bs.cell.mid,chia.2.bs.cell.mid,win.size=0.5e3,dists=NA,perms=100,...){
	if(i%%100==0) print(paste0('################### i: ',i))
	res <- lapply(1:perms,function(j) make_single_random_perm_ov_chia_pet_by_win(j,s.cell,p_chr=p_chr,e_chrs=e_chrs,
											prom.bs.mid,enh.bs.mid,chia.1.bs.cell.mid,chia.2.bs.cell.mid,dists=dists))
	res <- sapply(res,function(x) sum(x$eHit))
	return(res)
}

###############################################################################################################################################
# The function creates a distribution of distances b/w EP links per chromosome.
# The function retuns a vector of distances per chromosome and EP link to be used by matchIt.
###############################################################################################################################################
make_EP_dist_per_chr <- function(chr,g_to_e,prom.bs,enh.bs,...){
	print(paste0("Now working on chromosome: ", chr))
	g_list = names(g_to_e)
	g_list = strsplit(g_list,"[_:]")
	g_list = do.call(c,lapply(g_list,function(x) unlist(x)[2]))
	tmp_g_to_e <- g_to_e[which(g_list==chr)]
	ep_names <- do.call(c,lapply(names(tmp_g_to_e),function(x) paste0(x,'##',tmp_g_to_e[[x]])))
	
	shift.size=abs(start(enh.bs)-end(enh.bs))/2
	enh.bs.mid = shift(enh.bs,shift=shift.size)
	enh.bs.mid = promoters(enh.bs.mid,upstream=0,downstream=1)

	shift.size=abs(start(prom.bs)-end(prom.bs))/2
	prom.bs.mid = shift(prom.bs,shift=shift.size)
	prom.bs.mid = promoters(prom.bs.mid,upstream=0,downstream=1)
	
	ep_dist <- start(enh.bs.mid[unlist(tmp_g_to_e)])-rep(start(prom.bs.mid[names(tmp_g_to_e)]),times=sapply(tmp_g_to_e,length))
	names(ep_dist) <- ep_names
	return(ep_dist)
}

make_EP_dist_per_chr_fantom <- function(chr,g_to_e,prom.bs,enh.bs,...){
	print(paste0("Now working on chromosome: ", chr))
	g_list = names(g_to_e)
	g_list = strsplit(g_list,"\\$")
	g_list = do.call(c,lapply(g_list,function(x) unlist(x)[3]))
	tmp_g_to_e <- g_to_e[which(g_list==chr)]
	ep_names <- do.call(c,lapply(names(tmp_g_to_e),function(x) paste0(x,'##',tmp_g_to_e[[x]])))
	
	shift.size=abs(start(enh.bs)-end(enh.bs))/2
	enh.bs.mid = shift(enh.bs,shift=shift.size)
	enh.bs.mid = promoters(enh.bs.mid,upstream=0,downstream=1)

	shift.size=abs(start(prom.bs)-end(prom.bs))/2
	prom.bs.mid = shift(prom.bs,shift=shift.size)
	prom.bs.mid = promoters(prom.bs.mid,upstream=0,downstream=1)
	
	ep_dist <- start(enh.bs.mid[unlist(tmp_g_to_e)])-rep(start(prom.bs.mid[names(tmp_g_to_e)]),times=sapply(tmp_g_to_e,length))
	names(ep_dist) <- ep_names
	return(ep_dist)
}

###############################################################################################################################################
# The function creates a target vector of length(ep_links) with 1 marking a variable in target.gr and 0 o.w for MatchIt tool
###############################################################################################################################################
createMatchItTargetVector <- function(ep_dist,target_ep_links, ...){
	true.target <- match(target_ep_links,names(ep_dist))
	Tr <- rep(0,times=length(ep_dist))
	names(Tr) <- names(ep_dist)
	Tr[true.target] <- 1
	return(Tr)
}

###############################################################################################################################################
# The function creates the set of indices in Tr that constitutes the control (Crt) set using MatchIt tool
###############################################################################################################################################
createMatchItCrtSet <- function(ep_dist, Tr, ratio = 1, ...){
	Tr.new = Tr
	ep_dist.new = ep_dist
	if(length(Tr)>3e5){
		maxVal <- max(ep_dist[Tr==1])
		minVal <- min(ep_dist[Tr==1])
		id <- which(ep_dist<=maxVal & ep_dist>=minVal)
		Tr.new = Tr[id]
		ep_dist.new = ep_dist[id]
		
	}
	df <- as.data.frame(cbind(treatment=Tr.new,ep_dist=ep_dist.new))
	form <- paste0('treatment ~ ',paste(colnames(df)[-1],collapse=" + "))
	m.out <- matchit(as.formula(form), data = df, method = "nearest", ratio = ratio,replace=FALSE) 
	m.data <- match.data(m.out)
	#print(length(set.id.t <- rownames(m.data[m.data$treatment==1,])))
	set.id.c <- rownames(m.data[m.data$treatment==0,])
	#print(length(set.id.c))
	#print(length(set.id <- union(set.id.t,set.id.c)))
	return(set.id.c)
}

