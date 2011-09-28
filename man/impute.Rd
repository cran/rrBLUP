\name{impute}
\alias{impute}

\title{
Nearest-neighbor imputation
}
\description{
Imputes alleles based on the nearest-neighbor average. 
}
\usage{
impute(G, n.neighbor = "ALL", D = NULL)
}

\arguments{
  \item{G}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} markers.
}
  \item{n.neighbor}{
Number of nearest neighbors to use for imputation.  If not specified, the population mean is used.
}
  \item{D}{
Distance matrix (\eqn{n \times n}) between lines.  If not specified, calculated by \code{dist}.
}
}
\details{
Missing data in neighbors are omitted for the average.  If all of the neighbors are missing at any marker, an error is returned.  
This function assumes markers are coded numerically, e.g., \{aa,Aa,AA\} = \{-1,0,1\} for bi-allelic markers.  
Imputed alleles are fractional numbers.}
\value{
The imputed G matrix.
}
\examples{
#random population of 200 lines with 1000 markers
#5% missing data
G <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  for (j in 1:1000) {
    G[i,j] <- ifelse(runif(1)<0.05,NA,ifelse(runif(1)<0.5,-1,1))
  }
}
G.impute <- impute(G,n.neighbor=20)
}