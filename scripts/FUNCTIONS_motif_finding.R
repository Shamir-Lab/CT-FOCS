### Auxiliary functions for motif finding analysis ###
## count occurences of words in size wordsize for alphabet = {A,C,G,T}
## in sequence s
countOccurences <- function(i, s, wordsize=1, ...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	a <- s2c(s)
	seqinr::count(a, word = wordsize, alphabet = s2c("ACGT"))
}

## Create covariance matrix on the number of A,C,G,T and all pairwise letters in sequences 'seqs' for MatchIt tool
createMatchItCovMatrix <- function(seqs, ...){
	x.1 <- mclapply(1:length(seqs),function(i)  countOccurences(i,seqs[[i]]),mc.cores=10)
	x.2 <- mclapply(1:length(seqs),function(i)  countOccurences(i,seqs[[i]],wordsize=2),mc.cores=10)
	x.1 <- do.call(rbind,x.1)
	x.2 <- do.call(rbind,x.2)
	X <- cbind(x.1,x.2)
	rownames(X) <- 1:length(seqs)#names(seqs)
	rm(x.1);rm(x.2);gc()
	return(X)
}

## create a target vector of length(seqs) with 1 marking a variable in target.gr and 0 o.w for MatchIt tool
createMatchItTargetVector <- function(X,seqs,target.gr, ...){
	true.target <- match(names(target.gr),names(seqs))
	Tr <- rep(0,times=dim(X)[1])
	names(Tr) <- rownames(X)
	Tr[true.target] <- 1
	return(Tr)
}

## create the set of indices in Tr that constitutes the background (Bg) set using MatchIt tool
createMatchItBgSet <- function(X, Tr, ratio=1, distance = "mahalanobis", replace=FALSE, ...){
	df <- as.data.frame(cbind(treatment=Tr,X))
	form <- paste0('treatment ~ ',paste(colnames(df)[-1],collapse=" + "))
	err = rnorm(n = dim(df)[1], mean = 0, sd = 0.2)
	if(distance == "mahalanobis")
		df[,2:dim(df)[2]] <- df[,2:dim(df)[2]]+err
	m.out <- matchit(as.formula(form), data = df, method = "nearest", distance=distance, ratio = ratio,replace=replace) 
	m.data <- match.data(m.out)
	print(length(set.id.t <- as.numeric(rownames(m.data[m.data$treatment==1,]))))
	print(length(set.id.c <- as.numeric(rownames(m.data[m.data$treatment==0,]))))
	print(length(set.id <- union(set.id.t,set.id.c)))
	return(set.id.c)
}

############################## Known motif finding tools ####################################

## FIMO tool
## gr - GRanges object
## return GRanges of FIMO output
runFimo <- function(gr, pwm_db=NA, genome=genome, type = "encode", ...){
	seqs <- as.list(getSeq(genome, gr, as.character=TRUE))

	chr.names <- as.character(seqnames(gr))
	start.inds <- as.character(start(gr))
	end.inds <- as.character(end(gr))
	dgf.names <- paste(chr.names, start.inds, sep = ':')
	dgf.names <- paste(dgf.names, end.inds, sep = '-')
	
	write.fasta(sequences = seqs, names = dgf.names, file.out = paste0("../tmp/target_",type,".fasta"))


	fasta_get_markov_cmd <- 'fasta-get-markov '

	cmd_line <- paste0(fasta_get_markov_cmd,paste0("../tmp/target_",type,".fasta"), paste0(" ../tmp/bg_",type,".txt"))

	res <- system(cmd_line, intern = TRUE)

	fimo_cmd <- 'fimo'
	if(is.na(pwm_db))
		pwm_db <- '../meme/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme'
	out_folder <- paste0('--oc ../tmp/fimo_out_',type)
	bgfile <- paste0('--bgfile ../tmp/bg_',type,'.txt')
	fasfile <- paste0("../tmp/target_",type,".fasta")
	max_stored_scores <- '--max-stored-scores 100000'

	cmd_line <- paste(fimo_cmd,bgfile,out_folder,max_stored_scores,pwm_db,fasfile,collapse=" ")

	res <- system(cmd_line, intern = TRUE)

	# Import FIMO output

	fimo.gr <- import(con = paste0("../tmp/fimo_out_",type,"/fimo.gff"), format = "GFF")
	tokens <- strsplit(as.character(seqnames(fimo.gr)),"[:-]")

	seqnames <- sapply(tokens,function(x) x[1])
	starts <- sapply(tokens,function(x) as.numeric(x[2]))-1
	ends <- starts+end(fimo.gr)
	starts <- starts+start(fimo.gr)

	fimo.gr.fix <- GRanges(seqnames=seqnames, ranges=IRanges(start=starts,end=ends),strand=strand(fimo.gr))

	mcols(fimo.gr.fix) <- mcols(fimo.gr)
	fimo.gr.fix$orig_seqnames <- as.character(seqnames(fimo.gr))
	fimo.gr <- fimo.gr.fix
	rm(fimo.gr.fix);gc()

	fimo.gr$TF <- unlist(mclapply(fimo.gr$Name,function(x) unlist(strsplit(x,"_"))[1],mc.cores=10))

	return(fimo.gr)

}


findEnrichedTFsZtest <- function(fimo.gr, fimo.gr.crt, dgf.bs.reduced, dgf.bs.reduced.crt, tf.names, tf.widths, ...){
	p.values <- rep(1,length(tf.names))
	names(p.values) <- tf.names
	tg.hits <- rep(0,length(tf.names))
	names(tg.hits) <- tf.names
	tg.lens <- rep(0,length(tf.names))
	names(tg.lens) <- tf.names
	crt.hits <- rep(0,length(tf.names))
	names(crt.hits) <- tf.names
	crt.lens <- rep(0,length(tf.names))
	names(crt.lens) <- tf.names		
	total.width.tg <- sum(width(dgf.bs.reduced))
	total.width.crt <- sum(width(dgf.bs.reduced.crt))
	for (tf in tf.names){
		print(paste("handling tf:",tf))
		W <- tf.widths[tf]

		# if W>L_i remove L_i from L
		L.tg <- total.width.tg - sum(width(dgf.bs.reduced)[width(dgf.bs.reduced) < W])
		L.crt <- total.width.crt - sum(width(dgf.bs.reduced.crt)[width(dgf.bs.reduced.crt) < W])		
	
		x <- c(sum(tf == fimo.gr$TF),sum(tf == fimo.gr.crt$TF))
		n <- c(2*(L.tg - W + 1),2*(L.crt - W + 1))
		res <- prop.test(x, n, alternative = "greater")
		p.values[tf] <- res$p.value
		tg.hits[tf] <- x[1]
		tg.lens[tf] <- n[1]
		crt.hits[tf] <- x[2]
		crt.lens[tf] <- n[2]
	}
	return(list(pVal = p.values,tg.hits = tg.hits, tg.lens = tg.lens, crt.hits = crt.hits, crt.lens = crt.lens))
}

findEnrichedTFsHGtest <- function(fimo.gr, fimo.gr.crt, dgf.bs.reduced, dgf.bs.reduced.crt, tf.names, tf.widths, ...){
	p.values <- rep(1,length(tf.names))
	names(p.values) <- tf.names
	tg.hits <- rep(0,length(tf.names))
	names(tg.hits) <- tf.names
	tg.lens <- rep(0,length(tf.names))
	names(tg.lens) <- tf.names
	crt.hits <- rep(0,length(tf.names))
	names(crt.hits) <- tf.names
	crt.lens <- rep(0,length(tf.names))
	names(crt.lens) <- tf.names		
	total.width.tg <- sum(width(dgf.bs.reduced))
	total.width.crt <- sum(width(dgf.bs.reduced.crt))
	for (tf in tf.names){
		print(paste("handling tf:",tf))
		W <- tf.widths[tf]

		# if W>L_i remove L_i from L
		L.tg <- total.width.tg - sum(width(dgf.bs.reduced)[width(dgf.bs.reduced) < W])
		L.crt <- total.width.crt - sum(width(dgf.bs.reduced.crt)[width(dgf.bs.reduced.crt) < W])		
		x <- c(sum(tf == fimo.gr$TF),sum(tf == fimo.gr.crt$TF))
		y <- c(2*(L.tg - W + 1),2*(L.crt - W + 1))
		m <- sum(x) #number of white balls in the urn
		n <- sum(y)-m #number of black balls in the urn
		K <- y[1] #number of balls to be drawn
		q <- x[1] 
		tg.hits[tf] <- x[1]
		tg.lens[tf] <- y[1]
		crt.hits[tf] <- x[2]
		crt.lens[tf] <- y[2]
		res <- phyper(q-1, m, n, K, lower.tail=FALSE)
		p.values[tf] <- res
	}
	return(list(pVal = p.values,tg.hits = tg.hits, tg.lens = tg.lens, crt.hits = crt.hits, crt.lens = crt.lens))
}

computeHgPval <- function(dl,...){
	tf.names = dl$tf.names
	p.values <- rep(1,length(tf.names))
	names(p.values) <- tf.names
	for(tf in tf.names){
		x <- c(dl$tg.hits[tf],dl$crt.hits[tf])
		if(any(is.na(x))) next
		y <- c(2*dl$tg.lens[tf],2*dl$crt.lens[tf])
		m <- sum(x) #number of white balls in the urn
		n <- sum(y)-m #number of black balls in the urn
		K <- y[1] #number of balls to be drawn
		q <- x[1] 
		res <- phyper(q-1, m, n, K, lower.tail=FALSE)
		p.values[tf] <- res
	}
	return(p.values)
}

getNumOfPos <- function(fimo.gr, dgf.gr, tf.names, tf.widths, ...){
	tg.lens <- rep(0,length(tf.names))
	names(tg.lens) <- tf.names
	tg.hits <- rep(0,length(tf.names))
	names(tg.hits) <- tf.names
	total.width.tg <- sum(width(dgf.gr))
	for (tf in tf.names){
		print(paste("handling tf:",tf))
		W <- tf.widths[tf]
		# if W>L_i remove L_i from L
		y <- total.width.tg - sum(width(dgf.gr)[width(dgf.gr) < W])
		y <- 2*(y - W + 1)
		x <- sum(tf == fimo.gr$TF)
		tg.hits[tf] <- x
		tg.lens[tf] <- y
	}
	return(list(hits=tg.hits,lens=tg.lens,prop = tg.hits[tf.names]/tg.lens[tf.names]))
}


