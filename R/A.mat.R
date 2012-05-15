A.mat <- function(G,min.MAF=NULL,max.missing=NULL,impute=TRUE,tol=0.02,n.core=1,return.G=FALSE){
	
tcrossprod.na <- function(V,W) {
    crossprod.na <- function(x1, x2) {
        mean(x1 * x2, na.rm = TRUE)
    }
    n1 <- nrow(V)
    n2 <- nrow(W)
    result <- matrix(0,n1,n2)
    for (i in 1:n1) {
    	result[i,] <- apply(W,1,crossprod.na,V[i,])
    }
    return(result)
}
    
    
impute2 <- function(W, cov.mat, mean.vec) {
	n <- nrow(W)
	m <- ncol(W)
	S <- matrix(0,n,n)
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
			W[,i] <- Wi
		} else {D <- tcrossprod(Wi)}
		S <- S + D
	}	
	return(list(S=S,W.imp=W))
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

mono <- which(freq*(1-freq)==0)
G[,mono] <- 2*tcrossprod(one,matrix(freq[mono],length(mono),1))-1

freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
W <- (G[, markers] + 1 - 2 *freq.mat )/sig.A

if (!missing) {
	A <- tcrossprod(W)/m
	rownames(A) <- rownames(G)
	return(A)
} else {
    if (!impute) {
       if (n.core > 1) {
          library(multicore)
          it <- split(1:n,factor(cut(1:n,n.core,labels=FALSE)))
          A.list <- mclapply(it,function(ix){tcrossprod.na(W[ix,],W)})
          A <- numeric(0)
          for (i in 1:n.core) {A <- rbind(A,A.list[[i]])}
       } else {
          A <- tcrossprod.na(W,W)
	   }	
       rownames(A) <- rownames(G)
       return(A)
    } else {
    #impute = TRUE	
    isna <- which(is.na(W)) 
	W[isna] <- 0
	if (m < n) {
		print("There are fewer markers than lines: imputing with population means.")
	}
	if ((m < n)|(tol < 0)) {
		A <- tcrossprod(W)/m
		rownames(A) <- rownames(G)
		if (return.G) {
			G[,markers] <- W*sig.A - 1 + 2*freq.mat
			return(list(A=A,G.imp=G))		
		} else {
			return(A)
		}
	}

	#do EM algorithm
	mean.vec.new <- matrix(rowMeans(W),n,1)
	cov.mat.new <- cov(t(W))
    W[isna] <- NA
	A.new <- cov.mat.new + tcrossprod(mean.vec.new)
	err <- tol+1
	print("A.mat converging:")
	while (err >= tol) {
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
		W.imp <- numeric(0)
		for (i in 1:n.pieces) {
			S <- S + pieces[[i]]$S
			W.imp <- cbind(W.imp,pieces[[i]]$W.imp)
		}
		mean.vec.new <- matrix(rowMeans(W.imp),n,1)
		cov.mat.new <- (S-tcrossprod(mean.vec.new)*m)/(m-1)
		A.new <- cov.mat.new + tcrossprod(mean.vec.new)
		err <- norm(A.old-A.new,type="F")/n
		print(err,digits=3)
	}
	rownames(A.new) <- rownames(G)
	if (return.G) {
		G[,markers] <- W.imp*sig.A - 1 + 2*freq.mat
		return(list(A=A.new,G.imp=G))
	} else {
		return(A.new)
	}
  	} #else IMPUTE
} #else MISSING
} #A.mat

