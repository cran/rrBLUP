A.mat <- function(G,method="IBS",min.MAF=0.01,max.missing=1,PD=TRUE) {
	
	crossprod.na <- function(x1,x2) {
		mean(x1 * x2,na.rm=TRUE)		
	}
	tcrossprod.na <- function(G) {
	     apply(G,1,function(x,G){apply(G,1,crossprod.na,x)},G)
	}
	
	method <- toupper(method)
	n <- nrow(G)
	
	if (length(which(is.na(G))) > 0) {
	
	freq <- apply(G+1,2,function(x){mean(x,na.rm=TRUE)})/2
	frac.missing <- apply(G,2,function(x){length(which(is.na(x)))})/n
	markers <- which((freq >= min.MAF)&(1-freq >= min.MAF)&(frac.missing < max.missing))
  	freq <- freq[markers]
  	G <- G[,markers]
  	m <- ncol(G)
  
	if (method=="IBS") {
		A <- tcrossprod.na(G)+1
	} else if (method=="UAR") {
		X <- G + 1
		Y <- t(t(X)/sqrt(2*freq*(1-freq)))
		c <- sqrt(2*freq/(1-freq))
		z <- matrix(colMeans(t(Y)*c,na.rm=TRUE),n,1)
		one <- matrix(rep(1,n),n,1)
		q <- mean(c^2)
		A <- tcrossprod.na(Y)-tcrossprod(one,z)-tcrossprod(z,one)+tcrossprod(one)*q
	} else {stop("Invalid method")}

    if (PD) {
     as.matrix(nearPD(A)$mat)
	} else {A}

	} else {

	freq <- colMeans(G+1)/2
	markers <- which((freq >= min.MAF)&(1-freq >= min.MAF))
	freq <- freq[markers]
	G <- G[,markers]
	m <- ncol(G)
  
	if (method=="IBS") {
		A <- tcrossprod(G)/m+1
	} else if (method=="UAR") {
		X <- G + 1
		Y <- t(t(X)/sqrt(2*freq*(1-freq)))
		c <- sqrt(2*freq/(1-freq))
		z <- matrix(colSums(t(Y)*c),n,1)
		one <- matrix(rep(1,n),n,1)
		q <- sum(c^2)
		A <- (tcrossprod(Y)-tcrossprod(one,z)-tcrossprod(z,one)+tcrossprod(one)*q)/m
	} else {stop("Invalid method")}		

    A  
  } # if missing data	
}
