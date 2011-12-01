\name{mixed.solve}
\alias{mixed.solve}

\title{
Mixed-model solver
}
\description{
Calculates maximum-likelihood (ML/REML) solutions for mixed models of the form

\deqn{y = X \beta + Z u + \varepsilon}

where \eqn{\beta} is a vector of fixed effects and \eqn{u} is a vector of random effects with 
\eqn{Var[u] = K \sigma^2_u}.  The residual variance is \eqn{Var[\varepsilon] = I \sigma^2_e}.  This class
of mixed models, in which there is a single variance component other than the residual error,
has a close relationship with ridge regression (ridge parameter \eqn{\lambda = \sigma_e^2 / \sigma^2_u}).
}
\usage{
mixed.solve(y, Z=NULL, K=NULL, X=NULL, method="REML", 
        bounds=c(1e-09, 1e+09), SE=FALSE, return.Hinv=FALSE)
}

\arguments{
  \item{y}{
Vector (\eqn{n \times 1}) of observations.  Missing values (NA) are omitted, along with the corresponding rows of X and Z. 
}
  \item{Z}{
Design matrix (\eqn{n \times m}) for the random effects.  If not passed, assumed to be the identity matrix.
}
  \item{K}{
Covariance matrix (\eqn{m \times m}) for random effects; must be positive semi-definite.  If not passed, assumed to 
be the identity matrix.
}
  \item{X}{
Design matrix (\eqn{n \times p}) for the fixed effects.  If not passed, a vector of 1's is used
to model the intercept.
}
  \item{method}{
Specifies whether the full ("ML") or restricted ("REML") maximum-likelihood method is used.
}
  \item{bounds}{
Array with two elements specifying the lower and upper bound for the ridge parameter.
}
  \item{SE}{
If TRUE, standard errors are calculated.
}
  \item{return.Hinv}{
If TRUE, the function returns the inverse of \eqn{H = Z K Z' + \lambda I}.  This is useful for \code{\link{GWA}}.
}
}
\details{
This function can be used to predict marker effects or breeding values (see examples).  The numerical method 
is based on the spectral decomposition of \eqn{Z K Z'} and \eqn{S Z K Z' S}, where \eqn{S = I - X (X' X)^{-1} X'} is 
the projection operator for the nullspace of \eqn{X} (Kang et al., 2008). For marker effects where K = I, the function 
will run faster if K is not passed than if the user passes the identity matrix.
}
\value{
If SE=FALSE, the function returns a list containing
\describe{
\item{$Vu}{estimator for \eqn{\sigma^2_u}}
\item{$Ve}{estimator for \eqn{\sigma^2_e}}
\item{$beta}{estimator for \eqn{\beta}}
\item{$u}{BLUP(\eqn{u})}
\item{$LL}{maximized log-likelihood (full or restricted, depending on method)}
}

If SE=TRUE, the list also contains
\describe{
\item{$beta.SE}{standard error for \eqn{\beta}}
\item{$u.SE}{standard error for BLUP(\eqn{u}) \eqn{- u}}
}

If return.Hinv=TRUE, the list also contains
\describe{
\item{$Hinv}{the inverse of \eqn{H}}
}
}
\references{
Kang et al. 2008. Efficient control of population structure in model organism association mapping. 
Genetics 178:1709-1723.

Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255.
}

\examples{
#random population of 200 lines with 1000 markers
G <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  G[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

#random phenotypes
u <- rnorm(1000)
g <- as.vector(crossprod(t(G),u))
h2 <- 0.5  #heritability
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

#predict marker effects
ans <- mixed.solve(y,Z=G)  #By default K = I
accuracy <- cor(u,ans$u)

#predict breeding values
ans <- mixed.solve(y,K=A.mat(G))
accuracy <- cor(g,ans$u)

}

