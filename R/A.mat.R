A.mat <- function(G,min.MAF=0.01,n.adj=0,PD=TRUE) {
	
	crossprod.na <- function(x1,x2) {
		mean(x1 * x2,na.rm=TRUE)		
	}
	tcrossprod.na <- function(G) {
	    ncol(G)*apply(G,1,function(x,G){apply(G,1,crossprod.na,x)},G)
	}
	
    raw.UAR <- function(G,missing,min.MAF) {        
        n <- nrow(G)
        one <- matrix(rep(1,n),n,1)
         
        if (missing) {
            freq <- apply(G+1,2,function(x){mean(x,na.rm=TRUE)})/2
            markers <- which((freq >= min.MAF)&(1-freq >= min.MAF))
            m <- length(markers)
            X <- G[,markers] + 1 - 2*tcrossprod(one,matrix(freq[markers],m,1))
            var.A <- 2*sum(freq[markers]*(1-freq[markers]))
            return(tcrossprod.na(X)/var.A)
        } else {
            freq <- colMeans(G+1)/2
            markers <- which((freq >= min.MAF)&(1-freq >= min.MAF))
            m <- length(markers)
            X <- G[,markers] + 1 - 2*tcrossprod(one,matrix(freq[markers],m,1))
            var.A <- 2*sum(freq[markers]*(1-freq[markers]))
            return(tcrossprod(X)/var.A)
        }        
    } #raw.UAR
    
    b.yx <- function(y,x) {
    #formula for regression without intercept
        sum(y*x)/sum(x^2)
    }
    
	n <- nrow(G)
    m <- ncol(G)
    
	if (length(which(is.na(G))) > 0) {missing=TRUE} else {missing=FALSE}
	
    A.raw <- raw.UAR(G,missing=missing,min.MAF=min.MAF)
    
	if (n.adj == 0) {
        A <- A.raw
    } else {
        n.causal <- round(m/2)
        beta <- rep(0,n.adj)
        for (i in 1:n.adj) {
            ix <- sort(runif(m),index.return=TRUE)$ix
            causal.marker <- ix[1:n.causal]
            obs.marker <- ix[(n.causal+1):m]
            A.causal <- raw.UAR(G[,causal.marker],missing=missing,min.MAF=min.MAF)
            A.obs <- raw.UAR(G[,obs.marker],missing=missing,min.MAF=min.MAF)
            beta[i] <- b.yx(as.vector(A.causal-diag(n)),as.vector(A.obs-diag(n)))
        }
        A <- mean(beta)*A.raw + (1-mean(beta))*diag(n)
    }
    
    if (PD & missing) {return(as.matrix(nearPD(A)$mat))} else {return(A)}

}

