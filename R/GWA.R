GWA <-
function(y,M,Z=NULL,X = NULL,min.MAF=0.05,n.core=1,check.rank=FALSE) {
#assumes genotypes on [-1,1] scale
#missing data not allowed, impute first
#fractional genotypes OK

pi <- 3.14159
AS1 <- c(0.31938,-0.35656,1.78148,-1.82126,1.33027)  #for p-value approximation
AS2 <- 0.2316419

n <- length(y)
y <- matrix(y,n,1)
not.NA <- which(!is.na(y))

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
m <- ncol(M)  # number of markers
if (is.null(m)) {
  m <- 1
  M <- matrix(M,length(M),1)
}
n.line <- nrow(M)
if (is.null(Z)) {Z <- diag(n)}
stopifnot(nrow(Z)==n)
stopifnot(ncol(Z)==n.line)

Z <- as.matrix(Z[not.NA,])
X <- as.matrix(X[not.NA,])
n <- length(not.NA)
y <- matrix(y[not.NA],n,1)
Hinv <- mixed.solve(y,X=X,Z=Z,K=A.mat(M,n.core=n.core,min.MAF=min.MAF),return.Hinv=TRUE)$Hinv  
df <- p + 1

if (length(which(is.na(M))) > 0) {missing=TRUE} else {missing=FALSE}

score.calc <- function(M) {
  scores <- array(0,ncol(M))
  rownames(scores) <- colnames(M)
  for (i in 1:ncol(M)) {
    Mi <- M[,i]
  freq <- mean(Mi+1,na.rm=TRUE)/2
  MAF <- min(freq,1-freq)
  if (MAF < min.MAF) {
    scores[i] <- 0
  } else {

  if (missing) {
  NA.mark <- which(is.na(Mi))
  Mi[NA.mark] <- 0
  X2 <- cbind(X,Z%*%Mi)
  temp <- rep(0,n.line)
  temp[NA.mark] <- 1
  not.NA <- which(Z%*%temp!=1)
  n2 <- length(not.NA)
  X2 <- X2[not.NA,]
  y2 <- y[not.NA]
  H2inv <- Hinv[not.NA,not.NA]
  } else {
  n2 <- n
  X2 <- cbind(X,Z%*%Mi)
  y2 <- y
  H2inv <- Hinv
  }
  if (check.rank) {
    rXsnp <- qr(X2)$rank
  } else {
    rXsnp <- df
  }
  if (rXsnp != df) {
    scores[i] <- 0
  } else {
  W <- crossprod(X2,H2inv%*%X2)
  Winv <- solve(W)
  beta <- Winv %*% crossprod(X2,H2inv%*%y2)
  resid <- y2 - X2 %*% beta
  s2 <- as.double(crossprod(resid,H2inv%*%resid))/(n2-df)
  CovBeta <- s2*Winv
  F <- beta[df]^2/CovBeta[df,df]
  pvalue <- 1 - pf(F,1,n2-df)

  if (pvalue==0) {
  u <- 1/(1+AS2*sqrt(F))
  logp <- (-F/2-log(2*3.14159)/2+log(as.double(crossprod(AS1,c(u,u^2,u^3,u^4,u^5)))))/log(10)
  scores[i] <- -logp
  } else {  
  scores[i] <- -log(pvalue,10)
  } #if pvalue == 0

  } #ifelse Xsnp full rank
  } #if/else MAF < minMAF
  } #for i
    
  return(scores)
} #end score.calc

  if (n.core > 1) {
    it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
    library(multicore)
    scores <- unlist(mclapply(it,function(markers){score.calc(M[,markers])}))
   } else {
    scores <- score.calc(M)
   }      
  return(scores)
} #end function

