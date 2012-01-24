A.mat <- function (G, min.MAF = 0.01, n.core = 1, PD = TRUE) 
{
    crossprod.na <- function(x1, x2) {
        mean(x1 * x2, na.rm = TRUE)
    }
    tcrossprod.na <- function(Wk,W) {
        t(apply(Wk, 1, function(x,W) {apply(W, 1, crossprod.na, x)}, W))
    }

    if (length(which(is.na(G))) > 0) {missing=TRUE} else {missing=FALSE}
    n <- nrow(G)
    one <- matrix(1, n, 1)
    freq <- apply(G + 1, 2, function(x) {mean(x, na.rm = missing)})/2
    markers <- which((freq >= min.MAF) & (1 - freq >= min.MAF))
    m <- length(markers)
    W <- G[,markers] + 1 - 2 * tcrossprod(one, matrix(freq[markers],m, 1))
	if (n.core > 1) {
          library(multicore)
          it <- split(1:n,factor(cut(1:n,n.core,labels=FALSE)))
          A.list <- mclapply(it,function(ix){ if (missing) {tcrossprod.na(W[ix,],W)} else {tcrossprod(W[ix,],W)/m} })
          A <- numeric(0)
          for (i in 1:n.core) {A <- rbind(A,A.list[[i]])}
	} else {
	  if (missing) {A <- tcrossprod.na(W,W)} else {A <- tcrossprod(W)/m}
	}	
	var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
	A <- A/var.A
    if (PD & missing) {
    	if (length(which(is.na(A)))>0){
    		warning("A has missing values. Cannot coerce to PD.")
    		return(A)
    	} else {return(as.matrix(nearPD(A)$mat))}
    }
    else {
        return(A)
    }
}
