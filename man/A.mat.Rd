\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates the realized additive relationship matrix.
}
\usage{
A.mat(G,min.MAF=0.01,n.adj=0,PD=TRUE)
}

\arguments{
  \item{G}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\} = \{aa,Aa,AA\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency; default is 0.01.
}
\item{n.adj}{
Number of iterations for adjusting the matrix (see details).  Default is 0.
}
\item{PD}{
When there is missing data, the A matrix is not positive semidefinite by construction.  If PD=TRUE, the function \code{\link{nearPD}} is used to return the closest positive definite matrix, which can then be used in \code{\link{mixed.solve}}.
}
}
\details{
The A matrix is calculated as \eqn{W W'/c}, where \eqn{W_{ik} = G_{ik} + 1 - 2 p_k} and \eqn{p_k} is the frequency of the 1 allele at marker k.  The normalization constant is \eqn{c = 2 \sum_k {p_k (1-p_k)}}.  

The above formula is only exact as the number of markers tends to the size of the genome.  To correct for sampling error, an adjustment procedure similar to that proposed by Yang et al. (2010) can be used.  In each iteration, the marker set is randomly split in half, the A matrix is calculated for both subsets, and \eqn{vec(A1-I)} is regressed on \eqn{vec(A2-I)}.  If b denotes the average regression coefficient, the adjusted relationship matrix A* is calculated as \eqn{A^* = b(A-I) + I}. 

When alleles are missing, each pairwise calculation is based only on the markers that are present for both lines.  
}
\value{
\eqn{n \times n} additive relationship matrix
}
\references{
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
