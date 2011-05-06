\name{RR.BLUP.MR}
\alias{RR.BLUP.MR}

\title{
Ridge regression-BLUP with marker reduction
}
\description{
This function predicts complex traits via the mixed model

\deqn{y = X \beta + Z G u + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can include environmental factors and population structure.
The vector \eqn{u} represents a set of genetic markers, which are modeled as random effects with \eqn{Var[u] = K \sigma^2}. 
The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}.  This function combines ridge regression and 
subset selection, using restricted maximum-likelihood (\code{\link{RR.BLUP}}) to identify the ridge parameter 
(\eqn{\lambda = \sigma_e^2 / \sigma^2}), and cross-validation to identify the marker subset. The marker reduction 
procedure uses association analysis (\code{\link{GWA}}) to score and rank the markers, from which a series of models 
with progressively fewer markers is generated.
}
\usage{
RR.BLUP.MR(y, G, Z, X = NULL, K = NULL, Nfold = 3, 
           NumReduced = 8, MinMarkers = 10, SE = "FALSE")
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
Kinship matrix for the population; must be positive semidefinite.  If not passed, \eqn{K = G G' /m}.
}
  \item{Nfold}{
Number of folds to use in the cross-validation procedure
}
  \item{NumReduced}{
Number of reduced models to generate
}
  \item{MinMarkers}{
The series of reduced models is evenly spaced on a logarithmic scale from MinMarkers to the full marker set \eqn{m}.
}
  \item{SE}{
If "TRUE", the function will calculate standard errors.
}
}
\details{
For further details, see the reference by Endelman.
}
\value{
If SE="FALSE", the function returns a list containing
\describe{
\item{$beta}{the estimator for \eqn{\beta}}
\item{$u}{the BLUP solution for \eqn{u}.  Markers with identically zero effects were not present in
the optimal model.}
\item{$CV}{A data frame with the results of the cross-validation. The first column contains the number of markers
in each reduced model; the second column contains the average prediction accuracy across the folds.}
}

If SE="TRUE", the returned list also contains
\describe{
\item{$beta.SE}{the standard errors for the fixed effects}
\item{$u.SE}{the standard errors for the prediction error}
}
}
\references{
Endelman, J.B. (submitted) Coupling ridge regression-BLUP and association analysis to predict complex traits. 
}

\seealso{
\code{\link{RR.BLUP}}, \code{\link{GWA}}
}
\examples{
data(wheat.G)
data(wheat.y)
n <- nrow(wheat.G)
answer <- RR.BLUP.MR(wheat.y,G=wheat.G[,1:100],Z=diag(n),NumReduced=3)
}