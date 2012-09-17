\name{kinship.BLUP}
\alias{kinship.BLUP}

\title{
Genomic prediction by kinship-BLUP (deprecated)
}
\description{
***This function has been superseded by \code{\link{kin.blup}}; please refer to its help page.
}
\usage{
kinship.BLUP(y, G.train, G.pred=NULL, X=NULL, Z.train=NULL, 
     K.method="RR", n.profile=10, mixed.method="REML", n.core=1)
}

\arguments{
  \item{y}{
Vector (\eqn{n.obs \times 1}) of observations.  Missing values (NA) are omitted.
}
  \item{G.train}{
Matrix (\eqn{n.train \times m}) of unphased genotypes for the training population: \eqn{n.train} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\}; fractional (imputed) and missing (NA) alleles are allowed.  
}
  \item{G.pred}{
Matrix (\eqn{n.pred \times m}) of unphased genotypes for the prediction population: \eqn{n.pred} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\}; fractional (imputed) and missing (NA) alleles are allowed.  
}
  \item{X}{
Design matrix (\eqn{n.obs \times p}) of fixed effects.  If not passed, a vector of 1's is used 
to model the intercept.
}
  \item{Z.train}{
0-1 matrix (\eqn{n.obs \times n.train}) relating observations to lines in the training set.  If not passed
the identity matrix is used.
}
   \item{K.method}{
"RR" (default) is ridge regression, for which K is the realized additive relationship matrix computed with \code{\link{A.mat}}.  The option "GAUSS" is a Gaussian kernel (\eqn{K = e^{-D^2/\theta^2}}) and "EXP" is an exponential kernel (\eqn{K = e^{-D/\theta}}), where Euclidean distances \eqn{D} are computed with \code{\link{dist}}.
}

  \item{n.profile}{
For K.method = "GAUSS" or "EXP", the number of points to use in the log-likelihood profile for the scale parameter \eqn{\theta}.
}
  \item{mixed.method}{
Either "REML" (default) or "ML".
}
\item{n.core}{
For Mac, Linux, and UNIX users, setting n.core > 1 will enable parallel execution on a machine with multiple cores. R package multicore must be installed for this to work.  Do not run multicore from within the R GUI; you must use the command line.   
}
}
\value{
\describe{
\item{$g.train}{BLUP solution for the training set}
\item{$g.pred}{BLUP solution for the prediction set (when G.pred != NULL)}
\item{$beta}{ML estimate of fixed effects}
}
For GAUSS or EXP, function also returns
\describe{
\item{$profile}{log-likelihood profile for the scale parameter}
}
}
\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255.
}

\examples{
#random population of 200 lines with 1000 markers
G <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  G[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

#random phenotypes
g <- as.vector(crossprod(t(G),rnorm(1000)))
h2 <- 0.5 
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

#split in half for training and prediction
train <- 1:100
pred <- 101:200
ans <- kinship.BLUP(y=y[train],G.train=G[train,],G.pred=G[pred,],K.method="GAUSS")

#correlation accuracy
r.gy <- cor(ans$g.pred,y[pred])
}