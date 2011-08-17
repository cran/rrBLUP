mixed.solve <-
function(y,Z,K,X=NULL,method="REML",bounds=c(1e-9,1e9),SE=FALSE)  {

pi <- 3.14159

n <- length(y)  # number of observations
y <- matrix(y,n,1)
if (is.null(X)) {
  p <- 1
  X <- matrix(rep(1,n),n,1)
}
p <- ncol(X)  # number of fixed effects
if (is.null(p)) {
  p <- 1
  X <- matrix(X,length(X),1)
}
m <- ncol(Z)   # number of random effects
if (is.null(m)) {
  m <- 1
  Z <- matrix(Z,length(Z),1)
}

stopifnot(dim(Z)[1]==n)
stopifnot(dim(X)[1]==n)
stopifnot(dim(K)[1]==m)
stopifnot(dim(K)[2]==m)

Xt <- t(X)
Zt <- t(Z)
XtX <- crossprod(X,X)

r <- qr(XtX)$rank  # get rank(X'X) = rank(X) = rank(X') 
stopifnot(p==r)  #must be full rank

XtXinv <- solve(XtX)
A <- X %*% XtXinv %*% Xt
S <- diag(n) - A

offset <- sqrt(n)
Hb <- Z %*% K %*% Zt + offset*diag(n)  
Hb.system <- eigen(Hb,symmetric=TRUE)  
phi <- Hb.system$values - offset  #eigenvalues
U <- Hb.system$vectors 

SHbS <- S %*% Hb %*% S
SHbS.system <- eigen(SHbS,symmetric=TRUE)  
theta <- SHbS.system$values[1:(n-p)] - offset  #eigenvalues
Q <- SHbS.system$vectors[,1:(n-p)]

omega <- crossprod(Q,y)
omega.sq <- omega^2

if (method=="ML") {
  f.ML <- function(lambda,n,theta,omega.sq,phi) {
     n * log(sum(omega.sq/(theta+lambda))) + sum(log(phi+lambda)) 
  }
  soln <- optimize(f.ML,interval=bounds,n,theta,omega.sq,phi)
  lambda.opt <- soln$minimum
  df <- n
} else {
  # REML case
  f.REML <- function(lambda,n.p,theta,omega.sq) {
    n.p * log(sum(omega.sq/(theta+lambda))) + sum(log(theta+lambda))
  }
  soln  <- optimize(f.REML,interval=bounds,n-p,theta,omega.sq)
  lambda.opt <- soln$minimum
  df <- n-p
}

Vu.opt <- sum(omega.sq/(theta+lambda.opt))/df
Ve.opt <- lambda.opt * Vu.opt
H <- Hb + (lambda.opt - offset)*diag(n)
Hinv <- U %*% diag(1/(phi+lambda.opt)) %*% t(U)
beta <- solve(Xt%*%Hinv%*%X,Xt%*%Hinv%*%y)
u <- K %*% Zt %*% Hinv %*% (y - X %*% beta)
LL = -0.5*(soln$objective + df + df*log(2*pi/df))

if (SE==FALSE) {
  list(Vu = Vu.opt, Ve = Ve.opt, beta = as.vector(beta), u = as.vector(u), LL = LL)
} else {
  #also include standard errors for beta and u
  beta.var <- Vu.opt * solve(Xt %*% Hinv %*% X)
  beta.SE <- sqrt(diag(beta.var))
  C <- K %*% Zt
  u.var <- Vu.opt*(K - C%*%Hinv%*%t(C) + C%*%Hinv%*%X%*%beta.var%*%Xt%*%Hinv%*%t(C))
  u.SE <- sqrt(diag(u.var))
  list(Vu = Vu.opt, Ve = Ve.opt, beta = as.vector(beta), beta.SE = as.vector(beta.SE), u = as.vector(u), u.SE = as.vector(u.SE), LL = LL)
}

}  #end of function

