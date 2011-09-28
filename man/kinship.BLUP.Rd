\name{kinship.BLUP}
\alias{kinship.BLUP}

\title{
Genomic prediction by kinship-BLUP
}
\description{
This function uses \code{\link{mixed.solve}} for genomic prediction based on the kinship between lines. 
}
\usage{
kinship.BLUP(y, G.train, G.pred=NULL, X=NULL, Z.train=NULL, 
        D = NULL, K.method="RR", n.profile=10, mixed.method="REML")
}

\arguments{
  \item{y}{
Vector (\eqn{n.obs \times 1}) of observations
}
  \item{G.train}{
Matrix (\eqn{n.train \times m}) of unphased genotypes for the training population: \eqn{n.train} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\} = \{aa,Aa,AA\}; fractional (imputed) alleles are allowed.  
}
  \item{G.pred}{
Matrix (\eqn{n.pred \times m}) of unphased genotypes for the prediction population: \eqn{n.pred} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\} = \{aa,Aa,AA\}; fractional (imputed) alleles are allowed.  
}
  \item{X}{
Design matrix (\eqn{n.obs \times p}) of fixed effects.  If not passed, a vector of 1's is used 
to model the intercept.
}
  \item{Z.train}{
0-1 matrix (\eqn{n.obs \times n.train}) relating observations to lines in the training set.  If not passed
the identity matrix is used.
}
 \item{D}{
An \eqn{(n.train + n.test) \times (n.train + n.test)} matrix of genotypic distances, normalized to the interval
[0,1]. If not passed, generated using \code{dist}.
}
  \item{K.method}{
Method for the kinship matrix: "RR" (default) is ridge regression (\eqn{K = G G'/m}), "GAUSS" uses a 
Gaussian spatial model (\eqn{K = e^{-D^2/a^2}}), "EXP" uses an exponential spatial model (\eqn{K = e^{-D/a}}), 
and "MR" is a marker reduction procedure based on association analysis (\code{\link{GWA}}).
}
 
  \item{n.profile}{
For K.method = "GAUSS" or "EXP", the number of points to use in the log-likelihood profile for the scale parameter \eqn{a}.
For K.method = "MR", the number of reduced models in the cross-validation accuracy profile.  
}
  \item{mixed.method}{
Either "REML" (default) or "ML".
}
}
\details{
Missing alleles (NA) are not allowed when K.method = "RR" or "MR" (\code{\link{impute}} first).  

The efficiency of repeated calls to kinship.BLUP (e.g., in cross-validation) for methods "GAUSS" or "EXP" can be improved 
by making a single call to \code{dist} and then passing appropriate rearrangements of the distance matrix, e.g., 
D[c(train,pred),c(train,pred)].
}
\value{
\describe{
\item{$g.train}{BLUP solution for the training set}
\item{$g.pred}{BLUP solution for the prediction set (when G.pred != NULL)}
\item{$beta}{ML estimate of fixed effects}
}
When K.method = "GAUSS" or "EXP", function also returns
\describe{
\item{$profile}{log-likelihood profile for the scale parameter}
}
When K.method = "MR", function also returns
\describe{
\item{$profile}{cross-validation accuracy of the reduced models}
}
}
\references{
Endelman, J.B. (submitted) rrBLUP: An R package for ridge regression and kernel methods for genomic selection.
}

\seealso{
\code{\link{mixed.solve}}
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
pred <- 1:100
ans <- kinship.BLUP(y=y[train],G.train=G[train,],G.pred=G[pred,],K.method="GAUSS")

#correlation accuracy
r.gy <- cor(ans$g.pred,y[pred])
}