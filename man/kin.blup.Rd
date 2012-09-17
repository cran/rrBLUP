\name{kin.blup}
\alias{kin.blup}

\title{
Genotypic value prediction based on kinship
}
\description{
Genotypic value prediction by G-BLUP, where the genotypic covariance G can be additive or based on a Gaussian kernel.
}
\usage{
kin.blup(data,geno,pheno,GAUSS=FALSE,K=NULL,fixed=NULL,covariate=NULL,
         PEV=FALSE,n.core=1,theta.seq=NULL,reduce=TRUE)
}

\arguments{
  \item{data}{
Data frame with columns for the phenotype, the genotype identifier, and any environmental variables.
}
\item{geno}{
Character string for the name of the column in the data frame that contains the genotype identifier.
}
\item{pheno}{
Character string for the name of the column in the data frame that contains the phenotype.
}
\item{GAUSS}{
To model genetic covariance with a Gaussian kernel, set GAUSS=TRUE and pass the Euclidean distance for K (see below).
}
\item{K}{
There are three options for specifying kinship:  
(1) If K=NULL, genotypes are assumed to be independent \eqn{(G=I \: V_g)}. 
(2) For breeding value prediction, set GAUSS=FALSE and use an additive relationship matrix for K to create the model \eqn{(G=K \: V_g)}. 
(3) For the Gaussian kernel, set GAUSS=TRUE and pass the Euclidean distance for K to create the model \eqn{G_{ij}=e^{-(K_{ij}/\theta)^2} \: V_g}.  
}
\item{fixed}{
An array of strings containing the names of columns that should be included as (categorical) fixed effects in the mixed model.
}
\item{covariate}{
An array of strings containing the names of columns that should be included as covariates in the mixed model.
}
\item{PEV}{
When PEV=TRUE, the function returns the prediction error variance for the genotypic values (\eqn{PEV_i = Var[g^*_i-g_i]}). 
}
\item{n.core}{
Specifies the number of cores to use for parallel execution of the Gaussian kernel method.  Requires R package multicore and cannot be used with the R GUI.
}
\item{theta.seq}{
The scale parameter for the Gaussian kernel is set by maximizing the restricted log-likelihood over a grid of values. By default, the grid is constructed by dividing the interval (0,max(K)] into 10 points.  Passing a numeric array to this variable (theta.seq = "theta sequence") will specify a different set of grid points (e.g., for large problems you might want fewer than 10).  
}
\item{reduce}{
When reduce=TRUE and there are repeated observations, kin.blup reduces the dimension of the problem to equal the number of lines (see details).  To disable this feature, set reduce=FALSE.
}
}
\details{
This function is a wrapper for \code{\link{mixed.solve}} and thus solves mixed models of the form:
\deqn{y = X \beta + [Z \: 0] g + \varepsilon}
where \eqn{\beta} is a vector of fixed effects, \eqn{g} is a vector of random genotypic values with covariance 
\eqn{G = Var[g]}, and the residuals are i.i.d. (\eqn{Var[\varepsilon] = I \: V_e}).  The design matrix for the genetic values has been partitioned to illustrate that not all lines need phenotypes (i.e., for genomic selection).  Unlike \code{\link{mixed.solve}}, this function does not return estimates of the fixed effects, only the BLUP solution for the genotypic values.  It was designed to replace \code{\link{kinship.BLUP}} and to relieve the user of having to explicitly construct design matrices.  Variance components are estimated by REML and BLUP values are returned for every entry in K, regardless of whether it has been phenotyped.  The rownames of K must match the genotype labels in the data frame for phenotyped lines; missing phenotypes (NA) are simply omitted.  

Unlike its predecessor, this function does not handle marker data directly.  For breeding value prediction, the user must supply a relationship matrix, which can be calculated from markers with \code{\link{A.mat}}.  For Gaussian kernel predictions, pass the Euclidean distance matrix for K, which can be calculated with \code{\link{dist}}.  

In the terminology of mixed models, both the "fixed" and "covariate" variables are fixed effects (\eqn{\beta} in the above equation): the former are treated as factors with distinct levels while the latter are continuous with one coefficient per variable.  If no fixed effects of any kind are specified, an overall mean is included in the model.  

The prediction error variance (PEV) is the square of the SE of the BLUPs (see \code{\link{mixed.solve}}) and can be used to estimate the accuracy of BLUP predictions.

When the number of observations exceeds the number of lines, as happens with replicated trials or measurements, genomic predictions can be done by explicitly modeling the environmental variables in a one-step approach.  For unbalanced trials, this accounts for differing degrees of replication but is more computationally intensive than a two-step model.  When reduce=TRUE, the original mixed model is multiplied by \eqn{L = (Z'Z)^{-1/2}Z'} to generate a new mixed model with dimension equal to the number of lines, i.i.d. residuals, and the same variance parameters as in the original model.  (In general, only linear combinations of the original fixed effects parameters \eqn{\beta} are estimable in the reduced model--this is handled by the software.)  The reduced problem is similar to a two-step prediction with rotation of the adjusted means (Schulz-Streeck et al. 2012), but it is one-step in the sense that variance components and BLUPs are only estimated once. 
}
\value{
If GAUSS = FALSE and PEV = FALSE, a list containing
\describe{
\item{$Vg}{REML estimate of the genetic variance}
\item{$Ve}{REML estimate of the error variance}
\item{$g}{BLUP solution for the genetic values}
}
If PEV = TRUE, the list also includes
\describe{
\item{$PEV}{Prediction error variance for the genetic values}
}
If GAUSS = TRUE, the list also includes
\describe{
\item{$profile}{the log-likelihood profile for the scale parameter in the Gaussian kernel}
}
}

\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. doi: 10.3835/plantgenome2011.08.0024

Schulz-Streeck, T., J.O. Ogutu and H.-P. Piepho. 2012. Comparisons of single-stage and two-stage approaches to genomic selection. Theoretical and Applied Genetics. doi: 10.1007/s00122-012-1960-1
}

\examples{
#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  M[i,] <- ifelse(runif(1000)<0.5,-1,1)
}
rownames(M) <- 1:200
A <- A.mat(M)

#random phenotypes
u <- rnorm(1000)
g <- as.vector(crossprod(t(M),u))
h2 <- 0.5  #heritability
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

data <- data.frame(y=y,gid=1:200)

#predict breeding values
ans <- kin.blup(data=data,geno="gid",pheno="y",K=A)
accuracy <- cor(g,ans$g)

}
