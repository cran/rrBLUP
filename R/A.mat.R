A.mat <- function(G,method,min.MAF=0.01) {
  method = toupper(method)
  freq <- colMeans(G+1)/2
  markers <- which((freq >= min.MAF)&(1-freq >= min.MAF))
  freq <- freq[markers]
  G <- G[,markers]
  m <- ncol(G)
  n <- nrow(G)
  
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
}
