GWA <-
function(y,G,Z=NULL,X = NULL,K=NULL,min.MAF=0.05,check.rank=FALSE) {
#assumes genotypes on [-1,1] scale
#missing data not allowed, impute first
#fractional genotypes OK

pi <- 3.14159
AS1 <- c(0.31938,-0.35656,1.78148,-1.82126,1.33027)  #for p-value approximation
AS2 <- 0.2316419

n <- length(y)
y <- matrix(y,n,1)
if (is.null(X)) {
  p <- 1
  X <- matrix(rep(1,n),n,1)
}
p <- ncol(X)
if (is.null(p)) {
  p <- 1
  X <- matrix(X,length(X),1)
}
stopifnot(nrow(X)==n)
m <- ncol(G)  # number of markers
if (is.null(m)) {
  m <- 1
  G <- matrix(G,length(G),1)
}
t <- nrow(G)

if (is.null(Z)) {Z <- diag(n)}

stopifnot(nrow(Z)==n)
stopifnot(ncol(Z)==t)

if (is.null(K)) {
  K <- tcrossprod(G,G)/m
}

stopifnot(nrow(K)==ncol(K))
stopifnot(nrow(K)==t)

out <- mixed.solve(y,X=X,Z=Z,K=K)  
H <- out$Ve/out$Vu*diag(n)+tcrossprod(Z%*%K,Z)

Hinv <- solve(H)
df <- p + 1

scores <- rep(0,m)
freq <- colMeans(G+1)/2
for (i in 1:m) {
  MAF <- min(freq[i],1-freq[i])
  if (MAF < min.MAF) {
    scores[i] <- 0
  } else {
  
  Xsnp <- cbind(X,Z%*%G[,i])

  if (check.rank==TRUE) {
    rXsnp <- qr(Xsnp)$rank
  } else {
    rXsnp <- df
  }
  if (rXsnp != df) {
    scores[i] <- 0 
  } else {
  A <- crossprod(Xsnp,Hinv%*%Xsnp)
  Ainv <- solve(A)
  beta <- Ainv %*% crossprod(Xsnp,Hinv%*%y)
  resid <- y - Xsnp %*% beta
  s2 <- as.double(crossprod(resid,Hinv%*%resid))/(n-df)
  CovBeta <- s2*Ainv
  F <- beta[df]^2/CovBeta[df,df]
  pvalue <- 1 - pf(F,1,n-df)

  if (pvalue==0) {
  u <- 1/(1+AS2*sqrt(F))
  logp <- (-F/2-log(2*3.14159)/2+log(as.double(crossprod(AS1,c(u,u^2,u^3,u^4,u^5)))))/log(10)
  scores[i] <- -logp
  } else {  
  scores[i] <- -log(pvalue,10)
  } #end pvalue == 0

  } #end ifelse Xsnp full rank
  } #end if/else MAF < minMAF
} #end for

scores
} #end function

