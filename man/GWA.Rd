\name{GWA}
\alias{GWA}

\title{
Genome-wide association analysis
}
\description{
Performs genome-wide association analysis based on the mixed model

\deqn{y = X \beta + Z u + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can model both environmental
factors and population structure.  The variable \eqn{u} models the genetic background of each
line as a random effect with \eqn{Var[u] = K \sigma^2}.  The residual variance is 
\eqn{Var[\varepsilon] = I \sigma_e^2}.
}
\usage{
GWA(y, G, Z, X = NULL, K = NULL, min.MAF = 0.01, 
     check.rank = "FALSE")
}

\arguments{
  \item{y}{
Vector (\eqn{n \times 1}) of observations
}
  \item{G}{
Matrix (\eqn{t \times m}) of genotypes for \eqn{t} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\} = \{aa,Aa,AA\}.
}
  \item{Z}{
0-1 matrix (\eqn{n \times t}) relating observations to lines
}
  \item{X}{
Design matrix (\eqn{n \times p}) for the fixed effects.  If not passed, a vector of 1's is used 
to model the intercept.
}
  \item{K}{
Kinship matrix for the population; must be positive semidefinite.  If not passed, \eqn{K = G G' / m}.
}
  \item{min.MAF}{
Specifies the minimum minor allele frequency (MAF).  If a marker has a MAF less than min.MAF,
it is assigned a zero score.  
}
  \item{check.rank}{
If "TRUE", function will check the rank of the augmented design matrix for each marker. Markers for
which the design matrix is singular are assigned a zero score.
}
}
\details{
This function uses the iterative, generalized least-squares approach of Kang et al. (2010).  

The use of a minimum MAF is typically adequate to ensure the problem is well-posed.  However,
if an error message indicates the problem is singular, set check.rank to "TRUE".  This will slow 
down the algorithm but should fix the error.
}
\value{
Returns \eqn{m \times 1} vector of the marker scores, which equal \eqn{-log_{10}}(p-value)
}
\references{
Kang et al. 2010. Variance component model to account for sample structure in genome-wide association studies.
Nat. Genet. 42:348-354.

Endelman, J.B. (submitted) Coupling ridge regression-BLUP and association analysis to predict complex traits. 
}

\seealso{
\code{\link{RR.BLUP}}
}
\examples{
data(wheat.G)
data(wheat.y)
n <- nrow(wheat.G)
scores <- GWA(wheat.y,G=wheat.G[,1:100],Z=diag(n))
}
