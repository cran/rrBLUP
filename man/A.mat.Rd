\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates an additive relationship matrix.
}
\usage{
A.mat(G,method="IBS",min.MAF=0.01,max.missing=1,PD=FALSE,n.iter=50)
}

\arguments{
  \item{G}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\} = \{aa,Aa,AA\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{method}{
Must be one of "IBS", "raw.UAR", or "adj.UAR" (see details).
}
\item{min.MAF}{
Minimum minor allele frequency; default is 0.01.
}
\item{max.missing}{
Maximum proportion of entries that can be missing for each marker.  Markers above this limit are not used.
}
\item{PD}{
For both UAR models, as well as with the IBS option when alleles are missing, the A matrix is not positive semidefinite by construction.  If PD=TRUE, the function \code{\link{nearPD}} is used to return the closest positive definite matrix, which can then be used in \code{\link{mixed.solve}}
}
\item{n.iter}{
Number of iterations for calculating the adjustment factor in method "adj.UAR" (see details).
}
}
\details{
In the identity-by-state (IBS) method, A = 2*IBS = 1 + \eqn{G G' / m }.  The assumption that the coefficient of coancestry equals the IBS is consistent with coalescence theory, in which all alleles are derived from a common ancestor but at different times in the past.  The IBS method is suitable for genomic prediction but not for estimating heritability from variance components.  For that purpose it is better to use the unified additive relationship (UAR) model, which is based on gametic correlations (Yang et al. 2010; Powell et al. 2010).  The off-diagonal elements of A in the raw UAR model are  

\deqn{A_{ij} = (1/m) \sum_k {\frac{(X_{ik} - 2 p_k)(X_{jk} - 2 p_k)}{2 p_k (1-p_k)} }}

and the diagonal elements are given by

\deqn{A_{ii} = 1 + (1/m) \sum_k {\frac{X_{ik}^2 - (1 + 2 p_k) X_{ik} + 2 p_k^2}{2 p_k (1-p_k)}}}

where \eqn{X = G + 1} and \eqn{p_k} is the frequency of the 1 allele at marker k.  The adjusted UAR model corrects for sampling error by randomly dividing the markers into two groups, one designated as causal loci and the other as observed markers.  The raw UAR model is calculated for both groups, and the causal A is then regressed on the observed A.  This procedure is repeated n.iter times and averaged to get the shrinkage factor.

When alleles are missing, each pairwise calculation is based only on the markers that are present for both lines.  
}
\value{
\eqn{n \times n} additive relationship matrix
}
\references{
Powell, J.E., P.M. Visscher, and M.E. Goddard. 2010. Reconciling the analysis of IBD and IBS in complex trait studies. Nat. Genet. 11:800-805.

Yang et al. 2010. Common SNPs explain a large proportion of the heritability for human height. Nat. Genet. 42:565-569.
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
