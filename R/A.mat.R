A.mat <- function(G,min.MAF=NULL,max.missing=NULL,tol=0.01,n.core=1){
	
impute2 <- function(W, cov.mat, mean.vec) {
	n <- nrow(W)
	m <- ncol(W)
	S <- matrix(0,n,n)
	W.sum <- rep(0,n)
	for (i in 1:m) {
		Wi <- matrix(W[,i],n,1)
		missing <- which(is.na(Wi))
		if (length(missing) > 0) {
			not.NA <- setdiff(1:n,missing)
			Bt <- solve(cov.mat[not.NA,not.NA],cov.mat[not.NA,missing])
			Wi[missing] <- mean.vec[missing] + crossprod(Bt,Wi[not.NA]-mean.vec[not.NA])
			C <- cov.mat[missing,missing] - crossprod(cov.mat[not.NA,missing],Bt)
			D <- tcrossprod(Wi)
			D[missing,missing] <- D[missing,missing] + C
		} else {D <- tcrossprod(Wi)}
		W.sum <- W.sum + Wi
		S <- S + D
	}	
	return(list(S=S,W.sum=W.sum))
}

n <- nrow(G)
frac.missing <- apply(G,2,function(x){length(which(is.na(x)))/n})
missing <- max(frac.missing) > 0
freq <- apply(G + 1, 2, function(x) {mean(x, na.rm = missing)})/2
MAF <- apply(rbind(freq,1-freq),2,min)
if (is.null(min.MAF)) {min.MAF <- 1/(2*n)}
if (is.null(max.missing)) {max.missing <- 1 - 1/(2*n)}
markers <- which((MAF >= min.MAF)&(frac.missing <= max.missing)) 
m <- length(markers)
sig.A <- sqrt(2 * mean(freq[markers] * (1 - freq[markers])))
one <- matrix(1, n, 1)
freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
W <- (G[, markers] + 1 - 2 *freq.mat )/sig.A
if (!missing) {
	return(tcrossprod(W)/m)
} else {
    isna <- which(is.na(W)) 
	W[isna] <- 0
	if (m < n) {
		print("There are fewer markers than lines: imputing with population means.")
	}
	if ((m < n)|(tol < 0)) {return(tcrossprod(W)/m)}

	#do EM algorithm
	mean.vec.new <- matrix(rowMeans(W),n,1)
	cov.mat.new <- cov(t(W))
    W[isna] <- NA
	A.new <- cov.mat.new + tcrossprod(mean.vec.new)
	err <- tol+1
	count <- 0
	print("A.mat converging:")
	while (err >= tol) {
		count <- count + 1
		A.old <- A.new
		cov.mat.old <- cov.mat.new
		mean.vec.old <- mean.vec.new
		if (n.core > 1) {
            library(multicore)
			it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
			pieces <- mclapply(it,function(mark2){impute2(W[,mark2],cov.mat.old,mean.vec.old)})
		} else {
			pieces <- list()
			pieces[[1]] <- impute2(W,cov.mat.old,mean.vec.old)
		}
		n.pieces <- length(pieces)
		S <- matrix(0,n,n)
		W.sum <- rep(0,n)
		for (i in 1:n.pieces) {
			S <- S + pieces[[i]]$S
			W.sum <- W.sum + pieces[[i]]$W.sum
		}
		mean.vec.new <- matrix(W.sum/m,n,1)
		cov.mat.new <- (S-tcrossprod(mean.vec.new)*m)/(m-1)
		A.new <- cov.mat.new + tcrossprod(mean.vec.new)
		err <- norm(A.old-A.new,type="F")/n
		print(err,digits=3)
	}
	return(A.new)
}
}

