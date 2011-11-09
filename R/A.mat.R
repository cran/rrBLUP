A.mat <- function(G) {
  m <- ncol(G)
  n <- nrow(G)
  X <- G + 1
  freq <- colMeans(X/2)
  Y <- t(t(X)/sqrt(2*freq*(1-freq)))
  c <- sqrt(2*freq/(1-freq))
  z <- matrix(colSums(t(Y)*c),n,1)
  one <- matrix(rep(1,n),n,1)
  q <- sum(c^2)
  (tcrossprod(Y)-tcrossprod(one,z)-tcrossprod(z,one)+tcrossprod(one)*q)/m
}
