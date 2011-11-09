\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates an additive relationship matrix using the current population as the reference population.
}
\usage{
A.mat(G)
}

\arguments{
  \item{G}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\} = \{aa,Aa,AA\}.  Missing data are not allowed (\code{\link{impute}} first). 
}
}
\details{
This procedure calculates the additive relationship matrix as \eqn{A = W W'/m}, where \eqn{W_{ij} = (X_{ij} - 2 p_i)/\sqrt{2 p_i (1 - p_i)}} and \eqn{X = G + 1}.  This \eqn{A} matrix has the property that the average of its off-diagonal elements is zero, which can be interpreted as meaning the IBD relationships are relative to the current population as the reference population (see Powell et al. [2010] for more information).  Powell et al. (2010) suggest calculating the diagonal elements of \eqn{A} as 1 + F, where F is the correlation between the two gametes in a diploid individual.  This so-called raw UAR model, although unbiased, is in general not positive semi-definite.  The diagonal elements of \eqn{A = W W'/m} are slightly different than the raw UAR model, but the two models have the same average F for the population.
}
\value{
\eqn{n \times n} additive relationship matrix
}
\references{
Powell et al. 2010. Reconciling the analysis of IBD and IBS in complex trait studies. Nat. Genet. 11:800-805.
}
\examples{
#random population of 200 lines with 1000 markers
G <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  G[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

#Additive relationship matrix
A <- A.mat(G)

#Simulate phenotypes
Va <- 1 #additive genetic variance
g <- mvrnorm(mu=rep(0,200),Sigma=A*Va) 
h2 <- 0.5 #narrow-sense heritability
y <- g + rnorm(200,sd=sqrt((1-h2)/h2*Va))

#estimate heritability
ans <- mixed.solve(y=y,Z=diag(200),K=A)
h2.est <- ans$Vu/(ans$Vu+ans$Ve)

}
