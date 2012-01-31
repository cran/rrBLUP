\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates the realized additive relationship matrix.
}
\usage{
A.mat(G,min.MAF=NULL,max.missing=NULL,tol=0.01,n.core=1)
}

\arguments{
  \item{G}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\} = \{aa,Aa,AA\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency; default removes monomorphic markers.
}
\item{max.missing}{
Maximum proportion of missing data; default removes completely missing markers.
}
\item{tol}{
Specifies convergence criterion for EM algorithm when there is missing data (see details). If tol < 0, missing data are imputed with the population mean.
}
\item{n.core}{
For Mac, Linux, and UNIX users, setting n.core > 1 will enable parallel execution on a machine with multiple cores. R package multicore must be installed for this to work.  Do not run multicore from within the R GUI; you must use the command line.   
}
}
\details{
The A matrix is calculated as \eqn{W W'/c}, where \eqn{W_{ik} = G_{ik} + 1 - 2 p_k} and \eqn{p_k} is the frequency of the 1 allele at marker k.  The normalization constant is \eqn{c = 2 \sum_k {p_k (1-p_k)}}.  

When marker data are missing, an EM algorithm is used to impute genotypes and converge on the maximum likelihood solution for A.  The EM algorithm stops at iteration t when the RMS error = \eqn{n^{-1} \|A_{t} - A_{t-1}\|_2} < tol.  If the user passes a negative value for tol, or if the number of markers is less than the number of lines, missing alleles are imputed with the population mean for each marker.
}
\value{
\eqn{n \times n} additive relationship matrix
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
