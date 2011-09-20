mixed.solve <- function (y, Z, K, X = NULL, method = "REML", bounds = c(1e-09,1e+09), SE = FALSE) {
pi <- 3.14159
n <- length(y)
y <- matrix(y, n, 1)
if (is.null(X)) {
p <- 1
X <- matrix(rep(1, n), n, 1)
}
p <- ncol(X)
if (is.null(p)) {
p <- 1
X <- matrix(X, length(X), 1)
}
m <- ncol(Z)
if (is.null(m)) {
m <- 1
Z <- matrix(Z, length(Z), 1)
}
stopifnot(nrow(Z) == n)
stopifnot(nrow(X) == n)
stopifnot(nrow(K) == m)
stopifnot(ncol(K) == m)
Xt <- t(X)
Zt <- t(Z)
XtX <- crossprod(X, X)
rank.X <- qr(XtX)$rank
stopifnot(p == rank.X)
XtXinv <- solve(XtX)
S <- diag(n) - X %*% XtXinv %*% Xt
offset <- sqrt(n)
Hb <- Z %*% K %*% Zt + offset * diag(n)
Hb.system <- eigen(Hb, symmetric = TRUE)
phi <- Hb.system$values - offset
min.phi <- min(phi)
if (min.phi < 0) {warning(paste("K not positive semi-definite; min(phi) =",min.phi))}
U <- Hb.system$vectors
SHbS <- S %*% Hb %*% S
SHbS.system <- eigen(SHbS, symmetric = TRUE)
theta <- SHbS.system$values[1:(n - p)] - offset
Q <- SHbS.system$vectors[, 1:(n - p)]
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
W <- 0*diag(p)
if (p > 1) {
for (i in 1:(p-1)) {
for (j in (i+1):p) {
W[i,j] <- sum(UtX[,i]*UtX[,j]/(phi+lambda.opt))
W[j,i] <- W[i,j]
}
} #for i
}  #if p > 1
for (i in 1:p) {W[i,i] <- sum(UtX[,i]^2/(phi+lambda.opt))}
beta <- solve(W,crossprod(UtX,crossprod(U,y)/(phi+lambda.opt)))
u <- K %*% Zt %*% U %*% (crossprod(U,(y - X %*% beta))/(phi+lambda.opt))	    
LL = -0.5 * (soln$objective + df + df * log(2 * pi/df))
if (SE == FALSE) {
list(Vu = Vu.opt, Ve = Ve.opt, beta = as.vector(beta), u = as.vector(u), LL = LL)
} else {
beta.var <- Vu.opt * solve(W)
beta.SE <- sqrt(diag(beta.var))
C <- K %*% Zt
CU <- C %*% U
WW <- rep(0,m)
for (i in 1:m) {WW[i]<-sum(CU[i,]^2/(phi+lambda.opt))}
WWW <- matrix(rep(0,m*p),m,p)
for (i in 1:m) {
for (j in 1:p) {
WWW[i,j] <- sum(CU[i,]*UtX[,j]/(phi+lambda.opt))
}
}  #for i 
u.SE <- sqrt(Vu.opt * (diag(K) - WW + diag(WWW %*% beta.var %*% t(WWW))))
list(Vu = Vu.opt, Ve = Ve.opt, beta = as.vector(beta), beta.SE = as.vector(beta.SE), u = as.vector(u), u.SE = as.vector(u.SE), LL = LL)
}
}