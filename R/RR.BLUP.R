RR.BLUP <-
function(y,X,Z,K,method="REML",bounds=c(1e-9,1e9),SE="FALSE")  {

pi <- 3.14159

n <- length(y)  # number of observations
p <- dim(X)[2]  # number of fixed effects
if (is.null(p)) {
  p <- 1
  X <- matrix(X,length(X),1)
}
m <- dim(Z)[2]  # number of random effects
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
Hb.system <- svd(Hb)  
phi <- Hb.system$d - offset  #eigenvalues
U <- Hb.system$u 

SHbS <- S %*% Hb %*% S
SHbS.system <- svd(SHbS)  
theta <- SHbS.system$d[1:(n-p)] - offset  #eigenvalues
Q <- SHbS.system$u[,1:(n-p)]

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

Vg.opt <- sum(omega.sq/(theta+lambda.opt))/df
Ve.opt <- lambda.opt * Vg.opt
H <- Hb + (lambda.opt - offset)*diag(n)
Hinv <- U %*% diag(1/(phi+lambda.opt)) %*% t(U)
beta <- solve(Xt%*%Hinv%*%X,Xt%*%Hinv%*%y)
u <- K %*% Zt %*% Hinv %*% (y - X %*% beta)
LL = -0.5*(soln$objective + df + df*log(2*pi/df))

if (SE=="FALSE") {
  list(Vg = Vg.opt, Ve = Ve.opt, beta = beta, u = u, LL = LL)
} else {
  #also include standard errors for beta and u
  beta.var <- Vg.opt * solve(Xt %*% Hinv %*% X)
  beta.SE <- sqrt(diag(beta.var))
  C <- K %*% Zt
  u.var <- Vg.opt*(K - C%*%Hinv%*%t(C) + C%*%Hinv%*%X%*%beta.var%*%Xt%*%Hinv%*%t(C))
  u.SE <- sqrt(diag(u.var))
  list(Vg = Vg.opt, Ve = Ve.opt, beta = beta, beta.SE = beta.SE, u = u, u.SE = u.SE, LL = LL)
}

}  #end of function

