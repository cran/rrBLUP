kinship.BLUP <- function(y,G.train,G.pred=NULL,X=NULL,Z.train=NULL,K.method="RR",n.profile=10,mixed.method="REML") {
#assumes genotypes coded on [-1,1] scale
#continuous values OK

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
   K <- A.mat(G)
   soln <- mixed.solve(y=y,X=X,Z=Z,K=K,method=mixed.method)
   if (n.pred > 0) {
     list(g.train=soln$u[1:n.train],g.pred=soln$u[n.train+1:n.pred],beta=soln$beta)
   } else {
     list(g.train=soln$u[1:n.train],beta=soln$beta)
   }
} else {
  if ((K.method != "EXP")&(K.method != "GAUSS")) {stop("Invalid K.method")}
  # "exp" or "gauss"
  theta <- setdiff(seq(0,1,length.out=n.profile+1),0)
  LL <- rep(0,n.profile)
  soln <- list()

  D <- as.matrix(dist(G))/2/sqrt(m)
  
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
    list(profile=cbind(theta,LL),g.train=g.train,g.pred=g.pred,beta=soln[[max.LL]]$beta)
  } else {
    list(profile=cbind(theta,LL),g.train=g.train,beta=soln[[max.LL]]$beta)
  }

} #if K.method
} #function 
