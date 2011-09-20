\name{distance}
\alias{distance}

\title{
Genetic distance matrix
}
\description{
Calculates the Euclidean distance between lines in genotype space.
}
\usage{
distance(G)
}

\arguments{
  \item{G}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines with \eqn{m} bi-allelic markers.  
Genotypes should be coded as \{-1,0,1\} = \{aa,Aa,AA\}, although fractional (imputed) values are allowed.
}
}
\details{
Genetic distances are calculated according to 
\deqn{D_{ij} = \sqrt{(1/4m)\Sigma_{k}(G_{ik}-G_{jk})^2}}
Missing data (NA) are omitted.  

In a homozygous population, \eqn{1-D^2} is the identity-by-state kinship matrix.
}
\value{
An \eqn{n \times n} matrix of distances, normalized to the interval [0,1].
}
\references{
Piepho, H.P. 2009. Ridge regression and extensions for genomewide selection in maize. Crop Science 49:1165-1176.
}

\examples{
#random population of 200 lines with 1000 markers
G <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  G[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

D <- distance(G)
}
