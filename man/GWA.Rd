\name{GWA}
\alias{GWA}

\title{
Genome-wide association analysis
}
\description{
Performs genome-wide association analysis based on the mixed model

\deqn{y = X \beta + Z g + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can model both environmental
factors and population structure.  The variable \eqn{g} models the genetic background of each
line as a random effect with \eqn{Var[g] = K \sigma^2_g}.  The residual variance is 
\eqn{Var[\varepsilon] = I \sigma_e^2}.
}
\usage{
GWA(y, G, Z=NULL, X=NULL, K=NULL, min.MAF=0.05, 
     check.rank=FALSE)
}

\arguments{
  \item{y}{
Vector (\eqn{n \times 1}) of observations.  Missing values (NA) are omitted.
}
  \item{G}{
Matrix (\eqn{t \times m}) of genotypes for \eqn{t} lines with \eqn{m} bi-allelic markers.
Genotypes should be coded as \{-1,0,1\} = \{aa,Aa,AA\}.  Fractional (imputed) alleles are
allowed.  
}
  \item{Z}{
0-1 matrix (\eqn{n \times t}) relating observations to lines. If not passed, the identity matrix
is used.
}
  \item{X}{
Design matrix (\eqn{n \times p}) for the fixed effects.  If not passed, a vector of 1's is used 
to model the intercept.
}
  \item{K}{
Kinship matrix for the population; if not passed, the realized relationship matrix \eqn{G G'} is used.
}
  \item{min.MAF}{
Specifies the minimum minor allele frequency (MAF).  If a marker has a MAF less than min.MAF,
it is assigned a zero score.  
}
  \item{check.rank}{
If TRUE, function will check the rank of the augmented design matrix for each marker. Markers for
which the design matrix is singular are assigned a zero score.
}
}
\details{
This function implements the iterative, generalized least-squares method of Kang et al. (2010), using
\code{\link{mixed.solve}} for variance component estimation.

The use of a minimum MAF is typically adequate to ensure the problem is well-posed.  However,
if an error message indicates the problem is singular, set check.rank to TRUE.  This will slow 
down the algorithm but should fix the error.

Missing marker data are not allowed (\code{\link{impute}} first).
}
\value{
Returns \eqn{m \times 1} vector of the marker scores, which equal \eqn{-log_{10}}(p-value)
}
\references{
Kang et al. 2010. Variance component model to account for sample structure in genome-wide association studies.
Nat. Genet. 42:348-354.
}

\examples{
#random population of 200 lines with 1000 markers
G <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  G[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(t(G),u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

scores <- GWA(y=y,G=G)
}
