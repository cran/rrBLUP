kin.blup <- function(data,geno,pheno,GAUSS=FALSE,K=NULL,fixed=NULL,covariate=NULL,PEV=FALSE,n.core=1,theta.seq=NULL)  {

	make.full <- function(X) {
		svd.X <- svd(X)
		r <- max(which(svd.X$d>1e-8))
		return(as.matrix(svd.X$u[,1:r]))
	}
	
	names <- colnames(data)	
	ypos <- match(pheno,names)
	if (is.na(ypos)) {
		stop("Phenotype name does not appear in data.")
	} else {
		y <- as.numeric(data[,ypos])
	}

	not.miss <- which(!is.na(y))
	resid <- vector("numeric",length(y))*NA
	if (length(not.miss)<length(y)) {
		data <- data[not.miss,]
		y <- y[not.miss]
	} 
	n <- length(y)
   
	X <- matrix(1,n,1)		
	if (!is.null(fixed)) {
		p <- length(fixed)
		for (i in 1:p) {
			xpos <- match(fixed[i],names)
			xx <- factor(data[,xpos])	
			if (length(unique(xx)) > 1) {X <- cbind(X,model.matrix(~x-1,data.frame(x=xx)))}
		}
	}
	if (!is.null(covariate)) {
		p <- length(covariate)
		for (i in 1:p) {
			xpos <- match(covariate[i],names)
			X <- cbind(X,data[,xpos])
		}
	}	

	gid.pos <- match(geno,names)
	if (is.na(gid.pos)) {stop("Genotype name does not appear in data.")}

	not.miss.gid <- as.character(unique(data[,gid.pos]))
	if (is.null(K)) {
		gid <- not.miss.gid
		v <- length(gid)
		Z <- matrix(0,n,v)
		colnames(Z) <- gid
		Z[cbind(1:n,match(data[,gid.pos],gid))] <- 1
		X2 <- make.full(X)
		ans <- mixed.solve(y=y,X=X2,Z=Z,SE=PEV)
		resid[not.miss] <- y-X2%*%ans$beta-Z%*%ans$u
		if (PEV) {
			return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,PEV=ans$u.SE^2,resid=resid,pred=ans$u+as.numeric(colMeans(X2)%*%ans$beta)))
		} else {
			return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,resid=resid,pred=ans$u+as.numeric(colMeans(X2)%*%ans$beta)))
		}
	} else {

		if (inherits(K,what="dist")) {K <- as.matrix(K)}
		gid <- rownames(K)
		ix.pheno <- match(not.miss.gid,gid)
		miss.pheno.gid <- which(is.na(ix.pheno))
		if (length(miss.pheno.gid)>0) {
			stop(paste("The following lines have phenotypes but no genotypes:",paste(not.miss.gid[miss.pheno.gid],collapse=" ")))
		}
		miss.gid <- setdiff(gid,not.miss.gid)
		ix <- c(ix.pheno,match(miss.gid,gid))
		K <- K[ix,ix]
		v <- length(not.miss.gid)
		Z <- matrix(0,n,v)
		Z[cbind(1:n,match(data[,gid.pos],not.miss.gid))] <- 1
	
		X2 <- make.full(X)
		Z2 <- cbind(Z,matrix(0,n,nrow(K)-v))			
	
		if (!GAUSS) {
			ans <- mixed.solve(y=y,X=X2,Z=Z2,K=K,SE=PEV)
			ix <- match(gid,rownames(ans$u))
			resid[not.miss] <- y-X2%*%ans$beta-Z2%*%ans$u
			if (PEV) {
				return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid,pred=ans$u[ix]+as.numeric(colMeans(X2)%*%ans$beta)))
			} else {
				return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],resid=resid,pred=ans$u[ix]+as.numeric(colMeans(X2)%*%ans$beta)))
			}
		
		} else {

		if (is.null(theta.seq)) {
			theta <- setdiff(seq(0,max(K),length.out=11),0)
		} else {
			theta <- theta.seq
		}
		n.profile <- length(theta)		
		ms.fun <- function(theta) {
	    	soln <- list()
    		n.t <- length(theta)
	    	for (i in 1:n.t) {
   			   soln[[i]] <- mixed.solve(y=y,X=X2,Z=Z2,K=exp(-(K/theta[i])^2),SE=PEV)
    		}
	    	return(soln)
  		}

	if ((n.core > 1) & requireNamespace("parallel",quietly=TRUE)) {
    		it <- split(theta,factor(cut(theta,n.core,labels=FALSE)))
	    	soln <- unlist(parallel::mclapply(it,ms.fun,mc.cores=n.core),recursive=FALSE)
  		} else {
    		soln <- ms.fun(theta)
	  	}      

 		LL <- rep(0,n.profile)
	  	for (i in 1:n.profile) {LL[i] <- soln[[i]]$LL}
		ans <- soln[[which.max(LL)]]	
		profile <- cbind(theta,LL)
		ix <- match(gid,rownames(ans$u))
		resid[not.miss] <- y-X2%*%ans$beta-Z2%*%ans$u

		if (PEV) {
			return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid,pred=ans$u[ix]+as.numeric(colMeans(X2)%*%ans$beta)))
		} else {
			return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],resid=resid,pred=ans$u[ix]+as.numeric(colMeans(X2)%*%ans$beta)))
		}	

		} #else GAUSS
	} #else is.null(K)
} #kin.blup