\name{RR.BLUP}
\alias{RR.BLUP}

\title{
Ridge regression-BLUP
}
\description{
Calculates maximum-likelihood (ML/REML) solutions for mixed models of the form

\deqn{y = X \beta + Z u + \varepsilon}

where \eqn{\beta} is a vector of fixed effects and \eqn{u} is a vector of random effects with 
\eqn{Var[u] = K \sigma^2}.  The residual variance is \eqn{Var[\varepsilon] = I \sigma^2_e}.  This class
of mixed models, in which there are only two variance components (\eqn{\sigma^2} and \eqn{\sigma^2_e}), 
has a close relationship with ridge regression (shrinkage parameter \eqn{\lambda = \sigma_e^2 / \sigma^2}).
}
\usage{
RR.BLUP(y, X, Z, K, method = "REML", 
        bounds = c(1e-09, 1e+09), SE = "FALSE")
}

\arguments{
  \item{y}{
Vector (\eqn{n \times 1}) of observations
}
  \item{X}{
Design matrix (\eqn{n \times p}) for the fixed effects
}
  \item{Z}{
Design matrix (\eqn{n \times m}) for the random effects
}
  \item{K}{
Covariance matrix (\eqn{m \times m}) for random effects.  Must be positive semidefinite.
}
  \item{method}{
Specifies whether full ("ML") or restricted ("REML") maximum-likelihood method is used.
}
  \item{bounds}{
Array with two elements specifying the lower and upper bound for the ridge parameter.
}
  \item{SE}{
If "TRUE", the function will calculate standard errors.  
}
}
\details{
This function uses the spectral decomposition algorithm described by Kang et al. (2008).   
}
\value{
If SE="FALSE", the function returns a list containing
\describe{
\item{$Vg}{the estimator for \eqn{\sigma^2}}
\item{$Ve}{the estimator for \eqn{\sigma^2_e}}
\item{$beta}{the estimator for \eqn{\beta}}
\item{$u}{the BLUP solution for \eqn{u}}
\item{$LL}{the maximized log-likelihood (full or restricted, depending on method)}
}

If SE="TRUE", the returned list also contains
\describe{
\item{$beta.SE}{the standard errors for the fixed effects}
\item{$u.SE}{the standard errors for the prediction error}
}
}
\references{
Kang et al. 2008. Efficient control of population structure in model organism association mapping. 
Genetics 178:1709-1723.

Endelman, J.B. (submitted) Coupling ridge regression-BLUP and association analysis to predict complex traits. 
}

\examples{
data(wheat.G)
data(wheat.y)
n <- nrow(wheat.G)
answer <- RR.BLUP(wheat.y,X=rep(1,n),Z=wheat.G[,1:100],K=diag(100))
}

