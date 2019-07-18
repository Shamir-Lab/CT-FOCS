##### regression analysis functions ###############

fracIntronicEpLinks <- function(E_P,rg.cons,dreg.bs,...){
	id_genes <- match(names(E_P),names(rg.cons))
	id_rep <- rep(id_genes,times=sapply(E_P,length))
	rg.cons.ep <- rg.cons[id_rep]
	enh.ep <- dreg.bs[unlist(E_P)]
	
	res <- mclapply(seq_len(length(enh.ep)),function(i) countSingleIntronicEpLinks(i,enh.ep[i],rg.cons.ep[i]),mc.cores=20)
	return(length(which(res==T)))
}

countSingleIntronicEpLinks <- function(i,enh.bs,gene.bs,...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	hit.within <- findOverlaps(enh.bs,gene.bs,ignore.strand=T,select='first')
	return(!is.na(hit.within))
}

getEnhInWin <- function(i,tss.bs,enh.bs,win.size=1*10**6,....){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	chr <- as.character(seqnames(tss.bs))
	cand.enh <- which(seqnames(enh.bs)==chr)
	tss.bs.win <- promoters(tss.bs,upstream=win.size,downstream=win.size)
	hit.within <- findOverlaps(enh.bs[cand.enh],tss.bs.win,ignore.strand=T,select="first")
	cand.enh <- cand.enh[which(!is.na(hit.within))]
	tss.start <- start(promoters(tss.bs,upstream=0,downstream=1))
	dis <- abs(start(enh.bs[cand.enh])-tss.start)
	dis.ix <- sort(abs(dis), decreasing=F, index.return=T)$ix
	return(cand.enh[dis.ix])
}

getEnhInDom <- function(i,tss.bs,enh.bs,hic.bs.dom,....){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	chr <- as.character(seqnames(tss.bs))
	cand.enh <- which(seqnames(enh.bs)==chr)
	dom.ind <- findOverlaps(tss.bs,hic.bs.dom,ignore.strand=T,select="first")
	if(is.na(dom.ind))
		return(NA)
	hit.within <- findOverlaps(enh.bs[cand.enh],hic.bs.dom[dom.ind],ignore.strand=T,select="first")
	cand.enh <- cand.enh[which(!is.na(hit.within))]
	dis <- abs(start(enh.bs[cand.enh])-start(tss.bs))
	dis.ix <- sort(abs(dis), decreasing=F, index.return=T)$ix
	return(cand.enh[dis.ix])
}

fetch_regr_data<-function(ind,Mg,Me,g_to_e,k=10,...){
	y = Mg[ind,]
	es = g_to_e[[ind]]
#	if(length(es)>k)
#		es <- es[1:k]
	x = t(Me[es,])
	if(dim(x)[1]==1){
		x <- t(x)
		colnames(x) <- es
	}
	na_ids <- apply(x,2,function(col) any(is.na(col)))
	if(sum(na_ids)!=0) x <- x[,!na_ids]
	return(list(x=x,y=y))
}

################ Functions for running several models #####################
# The methods below return a list whose first three items are:
# 1. The regression model object
# 2. The p-value of the model vs. a model with the intercept only
# 3. A summary table of the coefficients

### Zero-inflated NB ####
run_zeroinf_test<-function(x,y,doEM=F,offs=NULL,...){
	if (all(y>0)){return (run_glmnb_test(x,y,offs=offs,...))}
	success = F
	try({
		if(is.null(offs)){
			zeroinf_res = zeroinfl(y ~ . | ., data = data.frame(x=x,y=y), dist = "negbin", EM = doEM,...)
		}
		else{
			zeroinf_res = zeroinfl(y ~ . | ., data = data.frame(x=x,y=y), dist = "negbin", EM = doEM,offset = offs,...)
		}
		success=T
	})
	if(!success ){
		return (list(reg_obj = NULL,pval = NA, coeffs = NULL,loglik = NA,dev=NA,success=F))
	}
	m0 = update(zeroinf_res,.~1)
	l1 = logLik(zeroinf_res);l0 = logLik(m0)
	df1 = attr(l1,"df");df0 = attr(l0,"df")
	pval = pchisq(2 * (l1[[1]] - l0[[1]]), df = df1-df0, lower.tail=FALSE)
	model_summary = summary(zeroinf_res)
	coeff_table = model_summary$coefficients$count
	return (list(reg_obj = zeroinf_res,pval = pval, coeffs = coeff_table,loglik = as.numeric(l1),dev=mean(zeroinf_res$residuals^2),success=T))
}
### Simple linear regression ####
run_lm_test<-function(x,y,...){
	lm_res = lm(y ~ . , data = data.frame(x=x,y=y))
	m0 = update(lm_res,.~1)
	l1 = logLik(lm_res);l0 = logLik(m0)
	df1 = attr(l1,"df");df0 = attr(l0,"df")
	logLik_pval = pchisq(2 * (l1[[1]] - l0[[1]]), df = df1-df0, lower.tail=FALSE)
	model_summary = summary(lm_res)
	fstats = model_summary$fstatistic
	pval = pf(fstats[1],fstats[2],fstats[3],lower.tail=F)
	coeff_table = model_summary[[4]]
	return (list(reg_obj = lm_res,pval = pval, 
		coeffs = coeff_table,logLik_pval=logLik_pval,r.squared = model_summary$r.squared,loglik = as.numeric(l1),dev=mean(lm_res$residuals^2),success = T))
}
### GLM NB regression ####
run_glmnb_test<-function(x,y,offs=NULL,...){
	success = F
	try({
		if(is.null(offs)){
			glm_res = glm.nb(y ~ ., data = data.frame(x=x,y=y),...)
		}
		else{
			glm_res = glm.nb(y ~ . + offset(offs), data = data.frame(x=x,y=y),...)
		}
		m0 = update(glm_res,.~1)
		l1 = logLik(glm_res);l0 = logLik(m0)
		df1 = attr(l1,"df");df0 = attr(l0,"df")
		pval = pchisq(2 * (l1[[1]] - l0[[1]]), df = df1-df0, lower.tail=FALSE)
		model_summary = summary(glm_res)
		coeff_table = model_summary$coefficients
		return (list(reg_obj = glm_res,pval = pval, coeffs = coeff_table,loglik = as.numeric(l1),dev=glm_res$deviance, success = T))
	})
	return (list(reg_obj = NULL,pval = NA, coeffs = NULL,loglik = NA,dev=NA,success=F))
}

### Computational permutation test ###
run_perm_test <- function(x,y,offs=NULL,func,repeats,statistic="loglik",...){
	fail_count <- 1
	res <- match.fun(func)(x=x,y=y,offs=offs)
	if(!res$success){
		return(NA)
	}
	stat_val_orig <- res[[statistic]]
	if(is.na(stat_val_orig)){
		return(NA)
	}
	for(i in seq_len(repeats)){
		try({
			curr_sample = sample(1:length(y))
			res <- match.fun(func)(x=x,y=y[curr_sample],offs=offs[curr_sample])
			stat_val <- res[[statistic]]
			if(is.na(stat_val)){
				next
			}
			if(statistic=="loglik"){
				if(stat_val>stat_val_orig)
					fail_count <- fail_count + 1
			}else
				if(stat_val<stat_val_orig)
					fail_count <- fail_count + 1
			}, TRUE)				
		
	}
	return(fail_count/(repeats+1))
}

################ Functions for stepwise validation #####################
run_single_batch_cv <- function(x,y,te_inds,func,offs=NULL,...){
	tr_x = x[-te_inds,];tr_y = y[-te_inds]
	te_x = x[te_inds,];te_y = y[te_inds]
	k <- dim(x)[2]
	if(k==1){
		tr_x <- t(t(tr_x))
		te_x <- t(t(te_x))
		colnames(tr_x) <- colnames(x)
		colnames(te_x) <- colnames(x)
	}
	te_preds <- rep(mean(tr_y),times=length(te_inds))#var(tr_y)==0 or model failed to run on training
	newdata = data.frame(te_x)
	try({
		if(is.null(dim(te_x)))
			newdata = data.frame(t(newdata))
		if(is.null(offs)){
			m  = func(tr_x,tr_y,...)
			colnames(newdata) <- tail(rownames(m$coeffs),k)
		}
		else{
			m  = func(tr_x,tr_y,offs=offs[-te_inds])
			colnames(newdata) <- rownames(m$coeffs)[2:(k+1)]
			newdata <- cbind(newdata,offs=offs[te_inds])
		}
		if(m$success)
			te_preds = predict(m[[1]],newdata = newdata ,type="response")
	})
	curr_samples = rownames(x)[te_inds]
	return(list(curr_samples=curr_samples,te_preds=te_preds))
}

run_batch_based_cross_validation_v2<-function(x,y,batches,func,offs=NULL,...){
	preds = c()
	all_batches = unique(batches)
	res <- mclapply(all_batches,function(b) run_single_batch_cv(x=x,y=y,te_inds=which(batches==b),func=func,offs=offs),mc.cores=10)
#	res <- sapply(all_batches,function(b) run_single_batch_cv(x=x,y=y,te_inds=which(batches==b),func=func,offs=offs))
	for(i in seq_len(length(res))){
		preds[res[[i]]$curr_samples] = res[[i]]$te_preds
	}
	preds = preds[rownames(x)]
	return(preds)
}

run_batch_based_cross_validation<-function(x,y,batches,func,offs=NULL,...){
	preds = c()
	all_batches = unique(batches)
	k <- dim(x)[2]
	for (b in all_batches){
		flag <- TRUE
		#if(b=="noCV") next;
		te_inds = which(batches==b)
		tr_x = x[-te_inds,];tr_y = y[-te_inds]
		te_x = x[te_inds,];te_y = y[te_inds]
		if(k==1){
			tr_x <- t(t(tr_x))
			te_x <- t(t(te_x))
			colnames(tr_x) <- colnames(x)
			colnames(te_x) <- colnames(x)
		}
		te_preds <- rep(mean(tr_y),times=length(te_inds))#var(tr_y)==0 or model failed to run on training
		newdata = data.frame(te_x)
		try({
			if(is.null(dim(te_x)))
				newdata = data.frame(t(newdata))
			if(is.null(offs)){
				m  = func(tr_x,tr_y,...)
				#surv_e <- unlist(strsplit(rownames(m$coeffs)[-1],"\\."))[seq(2,2*length(rownames(m$coeffs)[-1]),2)]
				suffix <- unlist(strsplit(rownames(m$coeffs)[-1],"\\."))[1] 
				colnames(newdata) <- paste0(suffix,".",colnames(newdata))
			}
			else{
				m  = func(tr_x,tr_y,offs=offs[-te_inds])
				suffix <- unlist(strsplit(rownames(m$coeffs)[-1],"\\."))[1] 
				colnames(newdata) <- paste0(suffix,".",colnames(newdata))
				newdata <- cbind(newdata,offs=offs[te_inds])
			}
			if(m$success)
				te_preds = predict(m[[1]],newdata = newdata ,type="response")
			flag <- FALSE
		})
		curr_samples = rownames(x)[te_inds]
		preds[curr_samples] = te_preds
		if(flag){
			print(b)
		}
	}
	preds = preds[rownames(x)]
	return(preds)
}

get_roc_eval<-function(preds,y_fpkm,threshold=1,...){
	if(sum(y_fpkm>threshold)==0){
		return (NA)		
	}
	if(sum(y_fpkm<threshold)==0){#all examples are positives - irrelevant to auc, because we need a mix of negatives/positives
		return (NA)		
	}
	curr_fpkm_factor = as.factor(as.numeric(y_fpkm>=threshold))
	return(auc(roc(preds,curr_fpkm_factor)))
}
get_spearman_corr<-function(preds,y_fpkm,threshold,return_pval=F,num.pos.th=3,...){
	if(sum(y_fpkm>threshold)==0){return (NA)}
	na_inds = is.na(preds) | is.na(y_fpkm)
	if(sum(na_inds) == length(preds)){return (NA)}
	preds = preds[!na_inds]; y_fpkm = y_fpkm[!na_inds]
	if(length(preds)<num.pos.th){return(NA)}
	curr_inds = y_fpkm>threshold
	if(return_pval){return(cor.test(preds[curr_inds],y_fpkm[curr_inds],method='spearman')$p.value)}
	return(cor(preds[curr_inds],y_fpkm[curr_inds],method='spearman'))
}
get_stepwise_validation_scores<-function(x,y,y_fpkm,threshold,batches,func,offs=NULL,...){
	preds = run_batch_based_cross_validation(x,y,batches,func,offs=offs,...)
	auroc = get_roc_eval(preds,y_fpkm,threshold,...)
	corr = get_spearman_corr(preds,y_fpkm,threshold,...)
	corr_pval = get_spearman_corr(preds,y_fpkm,threshold,T,...)
	return(list(auroc=auroc,corr=corr,corr_pval=corr_pval,preds=preds))
}

apply_cv_test_by_single_type <- function(j,arg_list,k=10,func = run_lm_test,pos.th=1,num.pos.th=3,...){
	if(j%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", j))}
	g_to_e <- arg_list[[5]]
	Mg_fpkm <- arg_list[[7]]
	p_name <- rownames(Mg_fpkm)[j]
	Mg <- arg_list[[1]]
	y_fpkm = Mg[p_name,]
	if(length(g_to_e[[p_name]]) < k || any(is.na(y_fpkm)) || sum(y_fpkm > pos.th) < num.pos.th){
			auroc = c(auroc =NA)
			corr = c(corr =NA)
			corr_pval = c(corr_pval=NA)
			preds = list(preds=rep(0,dim(Mg_fpkm)[2]))
			return(list(auroc,corr,corr_pval,preds))
	}
	Mg_normalized = arg_list[[2]]
	Me <- arg_list[[3]]
	Me_normalized <- arg_list[[4]]
	offs <- arg_list[[6]]
	col2type <- arg_list[[8]]
	#regr_data = fetch_regr_data(p_name,Mg,Me,g_to_e)
	regr_data_normalized = fetch_regr_data(p_name,Mg_normalized,Me_normalized,g_to_e)
	#x = regr_data$x;y=regr_data$y	
	x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
	if(ncol(x_n)<k){
			auroc = c(auroc =NA)
			corr = c(corr =NA)
			corr_pval = c(corr_pval=NA)
			preds = list(preds=rep(0,dim(Mg_fpkm)[2]))
			return(list(auroc,corr,corr_pval,preds))
	}
	if(ncol(x_n)>k) x_n <- x_n[,1:k]
	colnames(x_n) <- paste0('e',1:10)
	#colnames(x) <- paste0('e',1:10)
	flag=TRUE
	try({
	val1 = get_stepwise_validation_scores(x_n,y_n,y_fpkm,threshold=pos.th,col2type,func,num.pos.th,...)
	flag=FALSE
	})
	if(flag){
	  auroc = c(auroc =NA)
	  corr = c(corr =NA)
	  corr_pval = c(corr_pval=NA)
	  preds = list(preds=rep(0,dim(Mg_fpkm)[2]))
	  return(list(auroc,corr,corr_pval,preds))
	}
	auroc = c(auroc=val1[[1]])	
	corr = c(corr=val1[[2]])
	corr_pval = c(corr_pval=val1[[3]])
	preds = list(preds=val1[[4]])
	return(list(auroc,corr,corr_pval,preds))
}

apply_cv_test_by_type <- function(j,arg_list,k=10,...){
	if(j%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", j))}
	g_to_e <- arg_list[[5]]
	Mg_fpkm <- arg_list[[7]]
	y_fpkm = Mg_fpkm[j,]
		if(length(g_to_e[[j]]) < k || sum(y_fpkm > 1) < 3){
			auroc = c(zinb=NA,ols=NA,glmnb=NA)
			corr = c(zinb=NA,ols=NA,glmnb=NA)
			corr_pval = c(zinb=NA,ols=NA,glmnb=NA)
			preds = list(zinb=rep(0,dim(Mg_fpkm)[2]),ols=rep(0,dim(Mg_fpkm)[2]),glmnb=rep(0,dim(Mg_fpkm)[2]))
			return(list(auroc,corr,corr_pval,preds))
		}
	Mg <- arg_list[[1]]
	Mg_normalized = arg_list[[2]]
	Me <- arg_list[[3]]
	Me_normalized <- arg_list[[4]]
	offs <- arg_list[[6]]
	col2type <- arg_list[[8]]
	regr_data = fetch_regr_data(j,Mg,Me,g_to_e)
	regr_data_normalized = fetch_regr_data(j,Mg_normalized,Me_normalized,g_to_e)
	x = regr_data$x;y=regr_data$y	
	x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
	val1 = get_stepwise_validation_scores(x,y,y_fpkm,1,col2type,run_zeroinf_test,offs=offs )
	val2 = get_stepwise_validation_scores(x_n,y_n,y_fpkm,1,col2type,run_lm_test)
	val3 = get_stepwise_validation_scores(x,y,y_fpkm,1,col2type,run_glmnb_test,offs=offs)
	auroc = c(zinb=val1[[1]],ols=val2[[1]],glmnb=val3[[1]])	
	corr = c(zinb=val1[[2]],ols=val2[[2]],glmnb=val3[[2]])
	corr_pval = c(zinb=val1[[3]],ols=val2[[3]],glmnb=val3[[3]])
	preds = list(zinb=val1[[4]],ols=val2[[4]],glmnb=val3[[4]])
	return(list(auroc,corr,corr_pval,preds))
}

apply_cv_test_by_type_v2 <- function(j,g_to_e,y_fpkm,regr_data,regr_data_normalized,offs,col2type,k=10,...){
	if(j%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", j))}
	if(length(g_to_e) < k || sum(y_fpkm > 1) < 3){
			auroc = c(zinb=NA,ols=NA,glmnb=NA)
			corr = c(zinb=NA,ols=NA,glmnb=NA)
			corr_pval = c(zinb=NA,ols=NA,glmnb=NA)
			preds = list(zinb=rep(0,dim(Mg_fpkm)[2]),ols=rep(0,dim(Mg_fpkm)[2]),glmnb=rep(0,dim(Mg_fpkm)[2]))
			return(list(auroc,corr,corr_pval,preds))
	}
	x = regr_data$x;y=regr_data$y	
	x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
	val1 = get_stepwise_validation_scores(x,y,y_fpkm,1,col2type,run_zeroinf_test,offs=offs )
	val2 = get_stepwise_validation_scores(x_n,y_n,y_fpkm,1,col2type,run_lm_test)
	val3 = get_stepwise_validation_scores(x,y,y_fpkm,1,col2type,run_glmnb_test,offs=offs)
	auroc = c(zinb=val1[[1]],ols=val2[[1]],glmnb=val3[[1]])	
	corr = c(zinb=val1[[2]],ols=val2[[2]],glmnb=val3[[2]])
	corr_pval = c(zinb=val1[[3]],ols=val2[[3]],glmnb=val3[[3]])
	preds = list(zinb=val1[[4]],ols=val2[[4]],glmnb=val3[[4]])
	return(list(auroc,corr,corr_pval,preds))
}


### functions for pairwise correlations ###################################

get_intra_chr_enh <- function(enh.bs,rg.cons.gene,win.size=500000,...){
	rg.cons.gene.win <- trim(promoters(rg.cons.gene,upstream=win.size,downstream=win.size))
	hit.within <- findOverlaps(enh.bs,rg.cons.gene.win,ignore.strand = TRUE)
	match.id <- unique(queryHits(hit.within))
	return(match.id)
}

run_pairwise_corr <- function(y_fpkm,e_fpkm,...){
	return(cor(e_fpkm,y_fpkm,method='pearson'))
}

run_pairwise_corr_pv <- function(y_fpkm,e_fpkm,...){
	return(cor.test(e_fpkm,y_fpkm,method='pearson')$p.value)
}

run_pairwise_test_per_gene <- function(ind,y_fpkm,rg.cons.gene,Me_normalized,dreg.bs.mid,threshold=1,...){
	if(ind%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", ind))}
	tss_name <- rg.cons.gene$name
	chr <- seqnames(rg.cons.gene)
	enh_ids <- get_intra_chr_enh(dreg.bs.mid,rg.cons.gene,...)
	total_enh <- length(enh_ids)
	if(sum(y_fpkm>=threshold)==0 | total_enh==0){
		res <- list(num_enh=total_enh,enh_ids=enh_ids,cor_pval=rep(NA,total_enh),cor_val=rep(NA,total_enh))
		return(res)
	}
	cor_pval <- lapply(enh_ids,function(x) run_pairwise_corr_pv(y_fpkm,Me_normalized[x,]))
	cor_val <- lapply(enh_ids,function(x) run_pairwise_corr(y_fpkm,Me_normalized[x,]))
	return(list(num_enh=total_enh,enh_ids=enh_ids,cor_pval=unlist(cor_pval),cor_val=unlist(cor_val)))
}

correct_cor_pvals <- function(pairwise_cor,fdr_thr=1e-5,method='BH',...){
	vp <- unlist(sapply(pairwise_cor, function(x) x$cor_pval), use.names = F)
	vq = p.adjust(vp,method=method)
	thr = max(vp[vq<=fdr_thr],na.rm = T)
	n <- sum(unlist(sapply(pairwise_cor,function(x) x$num_enh)))
	pos=0
	for(i in seq_len(length(pairwise_cor))){
		corrected_pval <- vq[(pos+1):(pos+pairwise_cor[[i]]$num_enh)]
		passed_th_fdr <-  pairwise_cor[[i]]$cor_pval<=thr
		pairwise_cor[[i]]$corrected_pval <- corrected_pval
		pairwise_cor[[i]]$passed_th_fdr <- passed_th_fdr
		pos <- pos+pairwise_cor[[i]]$num_enh
		if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	}
	return(list(thr=thr,pairwise_cor=pairwise_cor))
}


get_sub_e_p_by_pearson_coef <- function(pairwise_cor,pe.r=0.7,...){
	res <- list()
	for(i in seq_len(length(pairwise_cor))){
		if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
		if(all(is.na(pairwise_cor[[i]]$passed_th_fdr)) | all(pairwise_cor[[i]]$passed_th_fdr==F))
			next
		enh_ids <- pairwise_cor[[i]]$enh_ids[pairwise_cor[[i]]$passed_th_fdr]
		above_th <- which(pairwise_cor[[i]]$cor_val[pairwise_cor[[i]]$passed_th_fdr]>=pe.r)
		if(length(above_th)==0)
			next
		enh_ids <- enh_ids[above_th]
		res[[names(pairwise_cor[i])]] <- pairwise_cor[[i]]
		res[[names(pairwise_cor[i])]]$enh_ids <- enh_ids
		res[[names(pairwise_cor[i])]]$num_enh <- length(enh_ids)
		res[[names(pairwise_cor[i])]]$model_ind <- i
	}
	return(res)

}

#Look at the OLS results: obtain estimation for the FDR and power of using R2
set_min_greater_than_zero<-function(x){
	if(sum(x==0,na.rm=T)==0){return(x)}
	x[x==0] = min(x[x>0],na.rm=T)
	return(x)
}

#correct enhancer p.vals by correction method - default "BY"
correct_p_matrix<-function(pmat,fdr_thr=0.05,method='BY'){
	vp = c(pmat)
	vq = p.adjust(vp,method=method)
	thr = max(vp[vq<=fdr_thr])
	return(list(thr=thr,corrected_mat = pmat<=thr))
}

#Get R.squared values
getRsquared <- function(i,y,y_pred,...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	return(cor(y,y_pred,method='pearson')**2)
}

#Get correlation values - either by Spearman or by Pearson (default)
getCor <- function(i,y,y_pred,...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	return(cor(y,y_pred,method='pearson'))
}

## proportional contribution r^2/R^2 (note that r is the Pearson coefficient)
get_pearson_coef <- function(i,ind,r2,Mg_normalized,Me_normalized,g_to_e,...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	regr_data_normalized = fetch_regr_data(ind,Mg_normalized,Me_normalized,g_to_e)
	x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
	return((cor(x_n,y_n,method='pearson')**2)/r2)
}

############################################################################################
############################## Linear mixed effect functions ###############################
############################################################################################
ind <- 1:k
form.glm <- paste0('y ~ ',paste('e',ind ,sep="",collapse=' + '))

calcEucDist <- function(exp.data){
  
  return(sort(rowSums(as.matrix(dist(t(exp.data),diag=T)/sqrt(dim(exp.data)[1])))/(dim(exp.data)[2]-1),decreasing=T))
  
}

fetch_regr_data<-function(ind,Mg,Me,g_to_e,k=10,...){
  y = Mg[ind,]
  es = g_to_e[[ind]]
  #	if(length(es)>k)
  #		es <- es[1:k]
  x = t(Me[es,])
  if(dim(x)[1]==1){
    x <- t(x)
    colnames(x) <- es
  }
  na_ids <- apply(x,2,function(col) any(is.na(col)))
  if(sum(na_ids)!=0) x <- x[,!na_ids]
  return(list(x=x,y=y))
}

# Chauvenet.R
#
# This function Chauvenet's criterion to remove outliers
Chauvenet <- function(datapoints, loop=TRUE, only.pos = FALSE, chau.pv.thr=0.5,...){
	numdatapoints <- nrow(datapoints)
	# calculate normalized distance to mean
	dist <- (datapoints - colMeans(datapoints))/apply(datapoints,2,sd)
	distAbs <- abs(dist)
	# calculate probability to have seen such a point assuming a normal
	# distribution
	prob <- apply(distAbs,c(1,2),function(x) numdatapoints*dnorm(x))
	# select only those points which have a probability >= 0.5
	sel <-  apply(prob,c(1,2),function(x) x>=chau.pv.thr)
	idx <- rowSums(sel) == ncol(datapoints)
	idx.out <- !idx
	neg.idx <- which(idx.out)
	if(only.pos & length(neg.idx)!=0){
	  neg.idx <- neg.idx[dist[neg.idx,1:dim(datapoints)[2], drop = FALSE]<0]
	  if(length(neg.idx)!=0){
	    idx.out[neg.idx] <- FALSE
	    idx[neg.idx] <- TRUE
	  }
	  
	}
	datapoints.out <- datapoints[idx.out,1:dim(datapoints)[2], drop = FALSE]
	datapoints <- datapoints[idx,1:dim(datapoints)[2], drop = FALSE]
	
	
	if(loop == TRUE){
		while(FALSE %in% idx){
      		numdatapoints <- nrow(datapoints)
      		dist <- (datapoints - colMeans(datapoints))/apply(datapoints,2,sd)
      		distAbs <- abs(dist)
      		prob <- apply(distAbs,c(1,2),function(x) numdatapoints*dnorm(x))
      		sel <-  apply(prob,c(1,2),function(x) x>=chau.pv.thr)
      		idx <- rowSums(sel) == ncol(datapoints)
      		if(!(FALSE %in% idx)) break
			    idx.out <- !idx
      		neg.idx <- which(idx.out)
      		if(only.pos & length(neg.idx)!=0){
      		  neg.idx <- neg.idx[dist[neg.idx,1:dim(datapoints)[2], drop = FALSE]<0]
      		  if(length(neg.idx)!=0){
      		    idx.out[neg.idx] <- FALSE
      		    idx[neg.idx] <- TRUE
      		  }
      		  
      		}
			    datapoints.out <- rbind(datapoints.out,datapoints[idx.out,1:dim(datapoints)[2], drop = FALSE])
      		datapoints <- datapoints[idx,1:dim(datapoints)[2], drop = FALSE]

    		}
	}
	
	if(dim(datapoints.out)[1]==0 & !only.pos){
		return(list('in'=datapoints,'out_pos'=datapoints.out,'out_neg'=datapoints.out))
	}
	if(dim(datapoints.out)[1]==0 & only.pos){
	  return(list('in'=datapoints,'out_pos'=datapoints.out,'out_neg'=datapoints[datapoints<0,,drop=F]))
	}
#	return(list('in'=datapoints,'out_pos'=datapoints.out[datapoints.out>0,,drop=F],'out_neg'=datapoints.out[datapoints.out<0,,drop=F]))
	if(only.pos)
	  return(list('out_pos'=datapoints.out,'out_neg'=datapoints[datapoints<0,,drop=F]))
	return(list('out_pos'=datapoints.out[datapoints.out>0,,drop=F],'out_neg'=datapoints.out[datapoints.out<0,,drop=F]))
}

outliersMAD <- function(data, MADCutOff = 3.5, bConstant = 1.4826, digits = 4, ...) {
   	#compute number of absolute MADs away for each value
    	#formula: abs( ( x - median(x) ) )/ mad(x)
	res <- matrix(0,nrow(data),ncol(data))
	rownames(res) <- rownames(data)
	colnames(res) <- colnames(data)
	mad_val <- apply(data,1,mad, constant = bConstant, na.rm = T)
	median_val <- apply(data,1,median,na.rm = T)
    	madAway <- t((t(data) - median_val)/mad_val)
    	#subset data that has absMADAway greater than the MADCutOff and replace them with replace
    	#can also replace values other than replace
	madAway[is.nan(madAway)] <- 0
    	res[madAway >= MADCutOff] <- 1
    	res[madAway <= MADCutOff*(-1)] <- (-1)
	return(list(madZ=round(madAway, digits),outliers=res))
}

getActiveP <- function(lmm_res, p_name, outliersMad, col2type, pv.thr=0.2, ...){
	is.sig=FALSE
	p.v.e <- lmm_res$p.v.e
	if(p.v.e[1]>pv.thr || !lmm_res$success) return(list(is.sig=is.sig,cells=NA))
	is.sig=TRUE
	cells <- unique(col2type[outliersMad[p_name,]==1])
	return(list(is.sig=is.sig,cells=cells))
}

getActiveEperP <- function(lmm_res, p_name, outliersMad, col2type, g_to_e, pv.thr=0.2, k=10, ...){
	p.v.e <- lmm_res$p.v.e
	sigE <- p.v.e[2:length(p.v.e)]<=pv.thr
	if(all(!sigE) || !lmm_res$success) return(list(cells=NA))
	eId <- g_to_e[[p_name]][1:k]
	names(eId) <- paste0('e',1:k)
	eId <- eId[sigE]
	eId <- eId[!is.na(eId)]
	cells <- lapply(eId, function(x) unique(col2type[outliersMad[x,]==1]))
	names(cells) <- names(eId)
	return(cells)
}

run_lme_test <- function(p_chr,dfs, k=10, method = "ML",...){
	success = F
	ind <- 1:k
	form.glm <- paste0('y ~ ',paste('e',ind ,sep="",collapse=' + '))
	try({
		GLM <- gls(as.formula(form.glm),data=dfs,method = method)
		lmm_res <- lme(as.formula(form.glm), random = ~ 1|cell, data = dfs, method = method, control = list(opt = "optim",optimMethod = "SANN"))
		res_glm=residuals(GLM)
		names(res_glm) <- dfs[names(res_glm),'cell']
		p.v.e.w = wilcox.test(abs(res_glm),abs(residuals(lmm_res)),paired=T,alternative="greater")$p.value	
		resid <- list(intersecpt=residuals(lmm_res))
		ranefs <- list(intersecpt=unlist(ranef(lmm_res)))
		p.v.e <- anova(GLM,lmm_res)$"p-value"[2]
		for(i in ind){
			lmm_res <- lme(as.formula(form.glm), random = list(cell=pdDiag(as.formula(paste0('~',paste('e',i,sep="",collapse='+'),' - 1')))), data = dfs, method = method, control = list(opt = "optim",optimMethod = "SANN"))
			p.v.e <- c(p.v.e,anova(GLM,lmm_res)$"p-value"[2])
			p.v.e.w = c(p.v.e.w,wilcox.test(abs(res_glm),abs(residuals(lmm_res)),paired=T,alternative="greater")$p.value)
			tmp <- list(residuals(lmm_res))
			names(tmp) <- paste0('e',i)
			resid <- c(resid,tmp)
			tmp <- list(unlist(ranef(lmm_res)))
			names(tmp) <- paste0('e',i)
			ranefs <- c(ranefs,tmp)
		}
		success=T
	})
	if(!success ){
		print(p_chr)
		return (list(p.v.e=rep(NA,k+1),resid=rep(list(null),k+1),ranefs=rep(list(null),k+1),success=F))
	}
	return(list(p.v.e=p.v.e,success=T,resid=resid,ranefs=ranefs))
}

run_glmer_nb_test <- function(p_chr,dfs.count, k=10,...){
	success = F
	ind <- 1:k
	form.glm <- paste0('y ~ ',paste('e',ind ,sep="",collapse=' + '))
	try({
		mad_res <- NULL
		glm_res <- glmmTMB(as.formula(paste0(form.glm, ' + offset(log(lib_size))')),
			data=dfs.count,
			ziformula=~1,
			family=nbinom2)
		glmm_res <- glmmTMB(as.formula(paste0(form.glm, ' + (1|cell)', ' + offset(log(lib_size))')),
			data=dfs.count,
			ziformula=~1,
			family=nbinom2)
		l1 = logLik(glmm_res);l0 = logLik(glm_res)
		df1 = attr(l1,"df");df0 = attr(l0,"df")
		pval = pchisq(2 * (l1[[1]] - l0[[1]]), df = df1-df0, lower.tail=FALSE)
		p.v.e <- pval
		mad_tmp <- outliersMAD(t(ranef(glmm_res)$cond$cell),...)
		mad_tmp <- colnames(mad_tmp$outliers)[mad_tmp$outliers==1]
		mad_res <- list(intercept=mad_tmp)
		for(i in ind){
			form.lmm.out <- paste0('y ~ ',paste('e',ind ,sep="",collapse=' + '), ' + (0 + ',paste('e',i,sep="",collapse='+'),'|cell)',' + offset(log(lib_size))')
			glmm_res <- glmmTMB(as.formula(form.lmm.out),
					data=dfs.count,
					ziformula=~1,
					family=nbinom2)
			l1 = logLik(glmm_res)
			df1 = attr(l1,"df")
			pval = pchisq(2 * (l1[[1]] - l0[[1]]), df = df1-df0, lower.tail=FALSE)
			p.v.e <- c(p.v.e,pval)
			mad_tmp <- outliersMAD(t(ranef(glmm_res)$cond$cell),...)
			mad_tmp <- colnames(mad_tmp$outliers)[mad_tmp$outliers==1]
			mad_l <- list(mad_tmp)
			names(mad_l) <- paste0('e',i) 
			mad_res <- c(mad_res,mad_l)
		}
		success=T
	})
	if(!success ){
		print(p_chr)
		return (list(p.v.e=rep(NA,k+1),mad_res=rep(list(NULL),k+1),success=F))
	}
	return(list(p.v.e=p.v.e,mad_res=mad_res,success=T))
}

#ind.sub should not include 0 (the random intercept), we always compute the random intercept only.
run_single_ranef_lmm_test <- function(p_chr,dfs, k=10, method="ML",...){
	success = F
	ind <- 1:k
	form.lmm.r <- paste0('~',paste('e',ind ,sep="",collapse='+'))
	form.glm <- paste0('y ~ ',paste('e',ind ,sep="",collapse=' + '))
	try({
		lmm_res <- NULL
		chau_res <- NULL
		lmm.null <- lme(as.formula(form.glm), random = ~ 1|cell, data = dfs, method = method, control = list(opt = "optim",optimMethod = "SANN"))
		lmm_res <- list(intercept=lmm.null)
		chau_tmp <- Chauvenet(ranef(lmm.null),...)
		chau_tmp <- rownames(chau_tmp$out_pos)
		chau_res <- list(intercept=chau_tmp)
		mad_tmp <- outliersMAD(t(ranef(lmm.null)),...)
		mad_tmp <- colnames(mad_tmp$outliers)[mad_tmp$outliers==1]
		mad_res <- list(intercept=mad_tmp)
		for(i in ind){
			lmm <- lme(as.formula(form.glm), random = list(cell=pdDiag(as.formula(paste0('~',paste('e',i,sep="",collapse='+'),' - 1')))), data = dfs, method = method, control = list(opt = "optim",optimMethod = "SANN"))
			lmm_l <- list(lmm)
			names(lmm_l) <- paste0('e',i) 
			lmm_res <- c(lmm_res,lmm_l)
			chau_tmp <- Chauvenet(ranef(lmm),...)
			chau_tmp <- rownames(chau_tmp$out_pos)
			chau_l <- list(chau_tmp)
			names(chau_l) <- paste0('e',i) 
			chau_res <- c(chau_res,chau_l)
			mad_tmp <- outliersMAD(t(ranef(lmm)),...)
			mad_tmp <- colnames(mad_tmp$outliers)[mad_tmp$outliers==1]
			mad_l <- list(mad_tmp)
			names(mad_l) <- paste0('e',i) 
			mad_res <- c(mad_res,mad_l)
		}
		success=T
	})
	if(!success ){
		print(p_chr)
		return (list(reg_obj = lmm_res,chau_res=chau_res,mad_res=mad_res,success=F))
	}
	return (list(reg_obj = lmm_res,chau_res=chau_res,mad_res=mad_res,success=T))
}


get_lmm_results <- function(i,p_chr,y_fpkm,Mg_normalized,Me_normalized,g_to_e,c2t,ranef.pv.thr=0.1,k=10,threshold=0,offs,...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
#  arguments <- list(...)
#  print(arguments)
	regr_data_normalized = fetch_regr_data(p_chr,Mg_normalized,Me_normalized,g_to_e)
	x_n = regr_data_normalized$x;y_n=unlist(regr_data_normalized$y)
	if(ncol(x_n)>k) x_n <- x_n[,1:k]
	cells <- unique(c2t)
	R <- matrix(0,nrow=dim(Mg_normalized)[2],ncol=length(cells))
	colnames(R) <- unique(cells)
	rownames(R) <- names(y_n)
	res <- t(sapply(seq_len(length(y_fpkm)),function(i){ if(y_fpkm[i]>threshold) R[i,c2t[i]]<-1; return(R[i,])}))
	rownames(res) <- names(y_n)
	e.names <- colnames(x_n)
	R <- res
	cell <- R
	x <- t(x_n)
	y <- y_n
	colnames(cell) <- make.names(colnames(cell))
	colnames(x) <- make.names(colnames(x))
	cell <- cell[,which(apply(cell,2,function(x) any(x==1)))]
	cell <- factor(c2t)
	names(cell) <-  make.names(c2t)
	x <- t(x)
	df = data.frame(x=x,y=y,cell=cell)
	names(df) <- c(paste0('e',1:k),'y','cell')
	dfs <- df
	dfs$lib_size=offs
	res <- run_lme_test(p_chr,dfs,k=k)
	p.v.e <- res$p.v.e
	resid <- res$resid
	ranefs <- res$ranefs
	res <- run_single_ranef_lmm_test(p_chr,dfs, k=k, method="ML",...)
	names(res$chau_res) <- c("intercept",e.names)
	names(res$mad_res) <- c("intercept",e.names)
	names(res$reg_obj) <- c("intercept",e.names)
	names(p.v.e) <- c("intercept",e.names)
	res[['p.v.e']] <- p.v.e
	res[['reg_obj']]<-NULL
	res[['resid']]<-resid
	res[['ranefs']]<-ranefs
	gc()
	return(res)
}


get_glmm_nb_results <- function(i,p_chr,y_fpkm,Mg_normalized,Me_normalized,g_to_e,c2t,ranef.pv.thr=0.1,k=10,threshold=0,offs,...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
#  arguments <- list(...)
#  print(arguments)
	regr_data_normalized = fetch_regr_data(p_chr,Mg_normalized,Me_normalized,g_to_e)
	x_n = regr_data_normalized$x;y_n=unlist(regr_data_normalized$y)
	if(ncol(x_n)>k) x_n <- x_n[,1:k]
	cells <- unique(c2t)
	R <- matrix(0,nrow=dim(Mg_normalized)[2],ncol=length(cells))
	colnames(R) <- unique(cells)
	rownames(R) <- names(y_n)
	res <- t(sapply(seq_len(length(y_fpkm)),function(i){ if(y_fpkm[i]>threshold) R[i,c2t[i]]<-1; return(R[i,])}))
	rownames(res) <- names(y_n)
	e.names <- colnames(x_n)
	R <- res
	cell <- R
	x <- t(x_n)
	y <- y_n
	colnames(cell) <- make.names(colnames(cell))
	colnames(x) <- make.names(colnames(x))
	cell <- cell[,which(apply(cell,2,function(x) any(x==1)))]
	cell <- factor(c2t)
	names(cell) <-  make.names(c2t)
	x <- t(x)
	df = data.frame(x=x,y=y,cell=cell)
	names(df) <- c(paste0('e',1:k),'y','cell')
	dfs <- df
	dfs$lib_size=offs
	res <- run_glmer_nb_test(p_chr,dfs,k=k,...)
	p.v.e <- res$p.v.e
	names(res$mad_res) <- c("intercept",e.names)
	names(p.v.e) <- c("intercept",e.names)
	res[['p.v.e']] <- p.v.e
	gc()
	return(res)
}

