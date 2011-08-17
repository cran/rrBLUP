kinship.BLUP <- function(y,G.train,G.pred=NULL,X=NULL,Z.train=NULL,K.method="RR",n.profile=10,mixed.method="REML") {
#assumes genotypes coded as {-1,0,1}

K.method <- toupper(K.method)

n.obs <- length(y)
y <- matrix(y,n.obs,1)

if (is.null(X)) {
  p <- 1
  X <- matrix(rep(1,n.obs),n.obs,1)
}
p <- ncol(X)  
if (is.null(p)) {
  p <- 1
  X <- matrix(X,length(X),1)
}

rX <- qr(X)$rank
stopifnot(rX==p)  #must be full rank design matrix
stopifnot(nrow(X)==n.obs)

if (is.null(Z.train)) {
  Z.train <- diag(n.obs)
}

m <- ncol(G.train)
n.train <- nrow(G.train)

stopifnot(ncol(Z.train)==n.train)
stopifnot(nrow(Z.train)==n.obs)

if (!is.null(G.pred)) {
  stopifnot(ncol(G.pred)==m)
  n.pred <- nrow(G.pred)
} else {
  n.pred <- 0
}

t <- n.pred + n.train #total number of lines

Z <- cbind(Z.train,matrix(rep(0,n.obs*n.pred),n.obs,n.pred))
G <- rbind(G.train,G.pred)

if (K.method == "RR") {
   freq <- colMeans(G+1)/2
   K.RR <- G %*% t(G)/2/sum(freq*(1-freq))
   soln <- mixed.solve(y=y,X=X,Z=Z,K=K.RR,method=mixed.method)
   if (n.pred > 0) {
     list(g.train=soln$u[1:n.train],g.pred=soln$u[n.train+1:n.pred],beta=soln$beta,Vg=soln$Vu,Ve=soln$Ve)
   } else {
     list(g.train=soln$u[1:n.train],beta=soln$beta,Vg=soln$Vu,Ve=soln$Ve)
   }
} else if (K.method == "MR") {
 
  #Partition lines in training set for cross-validation
  n.fold <- 5
  sets <- rep(0,n.train)
  ix <- sort(runif(n.train),index.return=TRUE)$ix
  SetSize <- floor(n.train/n.fold)
  for (k in 1:(n.fold-1)) {sets[ix[((k-1)*SetSize+1):(k*SetSize)]] <- k}
  sets[ix[((n.fold-1)*SetSize+1):n.train]] <- n.fold

  NumMarkers <- round(exp(seq(log(10),log(m),length.out=n.profile)))

  r.yy <- matrix(rep(0,n.fold*n.profile),n.profile,n.fold)
  for (j in 1:n.fold) {
    CV.test.lines <- which(sets==j)
    CV.train.lines <- setdiff(1:n.train,CV.test.lines)
    CV.test.obs <- NULL
    for (line in CV.test.lines) {
       CV.test.obs <- union(CV.test.obs,which(Z.train[,line]==1))
    }
    CV.train.obs <- setdiff(1:n.obs,CV.test.obs)
    CV.n.train <- length(CV.train.obs)
    CV.n.test <- length(CV.test.obs)

    Scores <- GWA(y=y[CV.train.obs],X=X[CV.train.obs,],Z=Z.train[CV.train.obs,CV.train.lines],G=G.train[CV.train.lines,])
    ScoreIdx <- sort(Scores,decreasing=TRUE,index.return=TRUE)$ix

    for (k in 1:n.profile) {
      MarkerList <- ScoreIdx[1:NumMarkers[k]]
      Zaug <- cbind(Z.train[CV.train.obs,CV.train.lines],matrix(rep(0,CV.n.train*CV.n.test),CV.n.train,CV.n.test))
      Gaug <- rbind(G.train[CV.train.lines,MarkerList],G.train[CV.test.lines,MarkerList])
      soln <- mixed.solve(y[CV.train.obs],X=X[CV.train.obs,],Z=Zaug,K=Gaug%*%t(Gaug),method=mixed.method)
      g.pred <- Z.train[CV.test.obs,CV.test.lines]%*%soln$u[CV.n.train+1:CV.n.test]
      r.yy[k,j] <- cor(g.pred,y[CV.test.obs])
    } #for k
  } #for j

  #find CV maximum
  r.yy.mean <- rowMeans(r.yy)
  opt.n.mark <- NumMarkers[which.max(r.yy.mean)]

  Scores <- GWA(y=y,X=X,Z=Z.train,G=G.train)
  ScoreIdx <- sort(Scores,decreasing=TRUE,index.return=TRUE)$ix
  MarkerList <- ScoreIdx[1:opt.n.mark]
  
  freq <- colMeans(G[,MarkerList]+1)/2
  K.RR <- G[,MarkerList] %*% t(G[,MarkerList])/2/sum(freq*(1-freq))
  soln <- mixed.solve(y=y,X=X,Z=Z,K=K.RR,method=mixed.method)

  if (n.pred > 0) {
    list(g.train=soln$u[1:n.train],g.pred=soln$u[n.train+1:n.pred],beta=soln$beta,Vg=soln$Vu,Ve=soln$Ve,profile=cbind(NumMarkers,r.yy.mean))
  } else {
    list(g.train=soln$u[1:n.train],beta=soln$beta,Vg=soln$Vu,Ve=soln$Ve,profile=cbind(NumMarkers,r.yy.mean))
  }

} else {
  # "exp" or "gauss"
  theta <- setdiff(seq(0,2,length.out=n.profile+1),0)
  LL <- rep(0,n.profile)
  soln <- list()

  D <- 0*diag(t)
  for (i in 1:(t-1)) {
   for (j in (i+1):t) {
      D[i,j] <- sqrt(sum((G[i,]-G[j,])^2))
   } #for j
  }  #for i

  D <- (D + t(D))/sqrt(m)
  for (i in 1:n.profile) {
    if (K.method == "EXP") {K <- exp(-D/theta[i])} 
    if (K.method == "GAUSS") {K <- exp(-(D/theta[i])^2) }
    soln[[i]] <- mixed.solve(y=y,X=X,Z=Z,K=K,method=mixed.method)
    LL[i] <- soln[[i]]$LL
  } #for i

  #find maximum LL soln
  max.LL <- which.max(LL)
  g.train <- soln[[max.LL]]$u[1:n.train]
  if (n.pred > 0) {
    g.pred <- soln[[max.LL]]$u[n.train+1:n.pred]
    list(profile=cbind(theta,LL),g.train=g.train,g.pred=g.pred,beta=soln[[max.LL]]$beta,Vg=soln[[max.LL]]$Vu,Ve=soln[[max.LL]]$Ve)
  } else {
    list(profile=cbind(theta,LL),g.train=g.train,beta=soln[[max.LL]]$beta,Vg=soln[[max.LL]]$Vu,Ve=soln[[max.LL]]$Ve)
  }

} #if K.method
} #function 

 


