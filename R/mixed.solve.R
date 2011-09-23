mixed.solve <- function (y, Z, K, X = NULL, method = "REML", bounds = c(1e-09,1e+09), SE = FALSE) {
pi <- 3.14159
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
m <- ncol(Z)
if (is.null(m)) {
m <- 1
Z <- matrix(Z,length(Z),1)
}
stopifnot(nrow(Z) == n)
stopifnot(nrow(X) == n)
stopifnot(nrow(K) == m)
stopifnot(ncol(K) == m)
XtX <- crossprod(X, X)
rank.X <- qr(XtX)$rank
stopifnot(p == rank.X)
XtXinv <- solve(XtX)
S <- diag(n) - tcrossprod(X%*%XtXinv,X)
if (n <= m) {
  spectral.method <- "eigen"
} else {
  spectral.method <- "cholesky"
  B <- try(chol(K),silent=TRUE)
  if (class(B)=="try-error") {
    stop("Error: K not positive definite.")
  }
}
if (spectral.method=="cholesky") {
ZBt <- tcrossprod(Z,B) 
svd.ZBt <- svd(ZBt,nu=n)
U <- svd.ZBt$u
phi <- c(svd.ZBt$d^2,rep(0,n-m))
SZBt <- S %*% ZBt
svd.SZBt <- svd(SZBt)
QR <- qr(cbind(X,svd.SZBt$u))
Q <- qr.Q(QR,complete=TRUE)[,(p+1):n]
R <- qr.R(QR)[(p+1):min(m+p,n),(p+1):(m+p)]
theta <- c(forwardsolve(t(R^2),svd.SZBt$d^2),rep(0,max(0,n-p-m)))
} else {
# spectral.method is "eigen"
offset <- sqrt(n)
Hb <- tcrossprod(Z%*%K,Z) + offset*diag(n)
Hb.system <- eigen(Hb, symmetric = TRUE)
phi <- Hb.system$values - offset
min.phi <- min(phi)
if (min.phi < -1e-6) {stop("Error: K is not positive semi-definite.")}
U <- Hb.system$vectors
SHbS <- S %*% Hb %*% S
SHbS.system <- eigen(SHbS, symmetric = TRUE)
theta <- SHbS.system$values[1:(n - p)] - offset
Q <- SHbS.system$vectors[, 1:(n - p)]
}  #if (n > m)
omega <- crossprod(Q, y)
omega.sq <- omega^2
if (method == "ML") {
f.ML <- function(lambda, n, theta, omega.sq, phi) {
 n * log(sum(omega.sq/(theta + lambda))) + sum(log(phi + lambda))
}
soln <- optimize(f.ML, interval = bounds, n, theta, omega.sq, phi)
lambda.opt <- soln$minimum
df <- n
} else {
f.REML <- function(lambda, n.p, theta, omega.sq) {
 n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
}
soln <- optimize(f.REML, interval = bounds, n - p, theta, omega.sq)
lambda.opt <- soln$minimum
df <- n - p
} #if method
Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
Ve.opt <- lambda.opt * Vu.opt
UtX <- crossprod(U,X)
W <- crossprod(UtX,UtX/(phi+lambda.opt))
beta <- solve(W,crossprod(UtX,crossprod(U,y)/(phi+lambda.opt)))
C <- K %*% crossprod(Z,U)
u <-  C %*% (crossprod(U,(y - X %*% beta))/(phi+lambda.opt))	    
LL = -0.5 * (soln$objective + df + df * log(2 * pi/df))
if (SE == FALSE) {
list(Vu = Vu.opt, Ve = Ve.opt, beta = as.vector(beta), u = as.vector(u), LL = LL)
} else {
beta.var <- Vu.opt * solve(W)
beta.SE <- sqrt(diag(beta.var))
WW <- tcrossprod(C%*%D,C)
WWW <- C%*%(UtX/(phi+lambda.opt))
u.SE <- sqrt(Vu.opt * (diag(K) - diag(WW) + diag(tcrossprod(WWW%*%beta.var,WWW))))
list(Vu = Vu.opt, Ve = Ve.opt, beta = as.vector(beta), beta.SE = as.vector(beta.SE), u = as.vector(u), u.SE = as.vector(u.SE), LL = LL)
}
}