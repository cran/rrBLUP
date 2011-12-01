\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates an additive relationship matrix.
}
\usage{
A.mat(G,method="IBS",min.MAF=0.01,max.missing=1,PD=TRUE)
}

\arguments{
  \item{G}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\} = \{aa,Aa,AA\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{method}{
Must be either "IBS" or "UAR", depending on what reference population is desired (see details).
}
\item{min.MAF}{
Minimum minor allele frequency; default is 0.01.
}
\item{max.missing}{
Maximum proportion of entries that can be missing for each marker.  Markers above this limit are not used.
}
\item{PD}{
This option only applies when alleles are missing, in which case A is not positive semidefinite by construction.  If PD=TRUE, the function \code{\link{nearPD}} is used to return the closest positive definite matrix, which can then be used in \code{\link{mixed.solve}}
}
}
\details{
For method "IBS", \eqn{A = 1 + G G' / m }, which means the coefficient of coancestry equals the identity-by-state (IBS).  The assumption that two alleles IBS are also identical by descent (IBD) is consistent with coalescence theory, in which all alleles are derived from a common ancestor but at different times in the past (Powell et al., 2010).

For method "UAR" (Unified Additive Relationship), \eqn{A = W W'/m}, where \eqn{W_{ij} = (X_{ij} - 2 p_i)/\sqrt{2 p_i (1 - p_i)}} and \eqn{X = G + 1}.  This \eqn{A} matrix has the property that the average of its off-diagonal elements is zero, which can be interpreted as meaning the IBD relationships are relative to the current population as the reference population (Powell et al., 2010).  In the raw UAR model of Powell et al. (2010), the diagonal elements of \eqn{A} are calculated as \eqn{1 + F}, where \eqn{F} is the correlation between the two gametes in a diploid individual.  This model, although unbiased, is in general not positive semi-definite.  The diagonal elements of \eqn{A = W W'/m} are slightly different than the raw UAR model of Powell et al. (2010), but the two models have the same average F for the population.

When alleles are missing, each pairwise calculation is based only on the markers that are present for both lines.
}
\value{
\eqn{n \times n} additive relationship matrix
}
\references{
Powell, J.E., P.M. Visscher, and M.E. Goddard. 2010. Reconciling the analysis of IBD and IBS in complex trait studies. Nat. Genet. 11:800-805.
}
\examples{
#random population of 200 lines with 1000 markers
G <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  G[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

#Additive relationship matrix
A <- A.mat(G)

}
