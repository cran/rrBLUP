\name{kinship.BLUP}
\alias{kinship.BLUP}

\title{
Genomic prediction by Kinship-BLUP
}
\description{
This function predicts complex traits via the mixed model

\deqn{y = X \beta + Z g + \varepsilon}

where \eqn{y} is a vector of observations for a training population, \eqn{X} is a full-rank design matrix for
the fixed effects \eqn{\beta}, \eqn{Z} is a design matrix relating lines to observations in the training set, and 
\eqn{g} is a vector of random effects for the genotypic values.  The covariance between lines is specified 
by \eqn{Var[g] = K \sigma^2_g}, with several options for the kinship matrix \eqn{K}.  
}
\usage{
kinship.BLUP(y, G.train, G.pred=NULL, X=NULL, Z.train=NULL, 
        K.method="RR", n.profile=10, mixed.method="REML")
}

\arguments{
  \item{y}{
Vector (\eqn{n.obs \times 1}) of observations
}
  \item{G.train}{
Matrix (\eqn{n.train \times m}) of genotypes for \eqn{n.train} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\} = \{aa,Aa,AA\}.
}
  \item{G.pred}{
Matrix (\eqn{n.pred \times m}) of genotypes for \eqn{n.pred} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\} = \{aa,Aa,AA\}.  
}
  \item{X}{
Design matrix (\eqn{n.obs \times p}) for the fixed effects.  If not passed, a vector of 1's is used 
to model the intercept.
}
  \item{Z.train}{
0-1 matrix (\eqn{n.obs \times n.train}) relating observations to lines in the training set.  If not passed
the identity matrix is used.
}
  \item{K.method}{
Method for the kinship matrix: "RR" (default) is the realized relationship matrix, "GAUSS" uses 
a Gaussian spatial model, "EXP" uses an exponential spatial model, and "MR" is a marker reduction procedure based
on association analysis (\code{\link{GWA}}).
}
  \item{n.profile}{
For K.method = "GAUSS" or "EXP", the number of points to use in the log-likelihood profile for the scale parameter.
For K.method = "MR", the number of reduced models in the cross-validation accuracy profile.  
}
  \item{mixed.method}{
Either "REML" (default) or "ML".
}
}
\details{
See the reference by Endelman.
}
\value{
\describe{
\item{$g.train}{BLUP solution for the training set}
\item{$g.pred}{BLUP solution for the prediction set (when G.pred != NULL)}
\item{$beta}{estimator for \eqn{\beta}}
\item{$Vg}{estimator for \eqn{\sigma^2_g}}
\item{$Ve}{estimator for \eqn{\sigma^2_e}}
}
When K.method = "GAUSS" or "EXP", predict() also returns
\describe{
\item{$profile}{log-likelihood profile for the scale parameter}
}
When K.method = "MR", predict() also returns
\describe{
\item{$profile}{cross-validation accuracy of the reduced models}
}
}
\references{
Endelman, J.B. (submitted) rrBLUP: An R package for ridge regression-BLUP and related genomic selection methods.
}

\seealso{
\code{\link{mixed.solve}}
}
\examples{
#random population of 500 inbred lines with 1000 markers
G <- matrix(rep(0,500*1000),500,1000)
for (i in 1:500) {
  G[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

#random phenotypes
g <- crossprod(t(G),rnorm(1000))
h2 <- 0.5 
y <- g + rnorm(500,mean=0,sd=sqrt((1-h2)/h2*var(g)))

#split in half for training and prediction
ans <- kinship.BLUP(y=y[1:250],G.train=G[1:250,],G.pred=G[251:500,])

#correlation accuracy
r.gy <- cor(ans$g.pred,y[251:500])
}