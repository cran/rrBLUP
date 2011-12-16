A.mat <- function(G,method="IBS",min.MAF=0.01,max.missing=1,PD=FALSE,n.iter=50) {
	
	crossprod.na <- function(x1,x2) {
		mean(x1 * x2,na.rm=TRUE)		
	}
	tcrossprod.na <- function(G) {
	     apply(G,1,function(x,G){apply(G,1,crossprod.na,x)},G)
	}
	
    to.vec <- function(A) {
        n <- nrow(A)
        const <- 100
        A.vec <- as.vector(triu(A,1)-const*tril(matrix(rep(1,n^2),n,n)))
        A.vec <- A.vec[-which(A.vec < (1 - const))]
        c(A.vec,diag(A)-1)
    }
    
    raw.UAR <- function(G,freq,missing,n) {        
        if (missing) {
            X <- G + 1
            Y <- t(t(X)/sqrt(2*freq*(1-freq)))
            c <- sqrt(2*freq/(1-freq))
            z <- matrix(colMeans(t(Y)*c,na.rm=TRUE),n,1)
            one <- matrix(rep(1,n),n,1)
            q <- mean(c^2)
            A <- tcrossprod.na(Y)-tcrossprod(one,z)-tcrossprod(z,one)+tcrossprod(one)*q
            for (i in 1:n) {A[i,i] <- 1 + mean((X[i,]^2-(1+2*freq)*X[i,]+2*freq^2)/(2*freq*(1-freq)),na.rm=TRUE)}         
        } else {
            X <- G + 1
            Y <- t(t(X)/sqrt(2*freq*(1-freq)))
            c <- sqrt(2*freq/(1-freq))
            z <- matrix(colSums(t(Y)*c),n,1)
            one <- matrix(rep(1,n),n,1)
            q <- sum(c^2)
            A <- (tcrossprod(Y)-tcrossprod(one,z)-tcrossprod(z,one)+tcrossprod(one)*q)/m
            for (i in 1:n) {A[i,i] <- 1 + mean((X[i,]^2-(1+2*freq)*X[i,]+2*freq^2)/(2*freq*(1-freq)))}
        }        
    A
    } #raw.UAR
    
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
	} else if (method=="RAW.UAR") {
        A <- raw.UAR(G,freq,missing=TRUE,n)
	} else if (method=="ADJ.UAR") {
        n.causal <- round(m/2)
        beta <- rep(0,n.iter)
        for (i in 1:n.iter) {
            ix <- sort(runif(m),index.return=TRUE)$ix
            causal.marker <- ix[1:n.causal]
            obs.marker <- ix[(n.causal+1):m]
            A.causal <- raw.UAR(G[,causal.marker],freq[causal.marker],missing=TRUE,n)
            A.obs <- raw.UAR(G[,obs.marker],freq[obs.marker],missing=TRUE,n)
            beta[i] <- cov(to.vec(A.causal),to.vec(A.obs))/var(to.vec(A.obs))     
        }
        A.raw <- raw.UAR(G,freq,missing=TRUE,n)
        A <- mean(beta)*A.raw
        A[cbind(1:n,1:n)] <- 1 + mean(beta)*(diag(A.raw)-1)
    } else {stop("Invalid method")}

    } else {
    # no missing data
	freq <- colMeans(G+1)/2
	markers <- which((freq >= min.MAF)&(1-freq >= min.MAF))
	freq <- freq[markers]
	G <- G[,markers]
	m <- ncol(G)
  
	if (method=="IBS") {
		A <- tcrossprod(G)/m+1
	} else if (method=="RAW.UAR") {
        A <- raw.UAR(G,freq,missing=FALSE,n)
    } else if (method=="ADJ.UAR") {
        n.causal <- round(m/2)
        beta <- rep(0,n.iter)
        for (i in 1:n.iter) {
            ix <- sort(runif(m),index.return=TRUE)$ix
            causal.marker <- ix[1:n.causal]
            obs.marker <- ix[(n.causal+1):m]
            A.causal <- raw.UAR(G[,causal.marker],freq[causal.marker],missing=FALSE,n)
            A.obs <- raw.UAR(G[,obs.marker],freq[obs.marker],missing=FALSE,n)
            beta[i] <- cov(to.vec(A.causal),to.vec(A.obs))/var(to.vec(A.obs))     
        }
        A.raw <- raw.UAR(G,freq,missing=FALSE,n)
        A <- mean(beta)*A.raw
        A[cbind(1:n,1:n)] <- 1 + mean(beta)*(diag(A.raw)-1)
    } else {stop("Invalid method")}		
  } # if missing data	

if (PD) {
     as.matrix(nearPD(A)$mat)
  } else {A}
}
