\name{GWAS}
\alias{GWAS}

\title{
Genome-wide association analysis
}
\description{
Performs genome-wide association analysis based on the mixed model (Yu et al. 2006):

\deqn{y = X \beta + Z g + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}.
}
\usage{
GWAS(pheno, geno, fixed=NULL, K=NULL, n.PC=0,
     min.MAF=0.05, n.core=1, P3D=TRUE, plot=TRUE)
}

\arguments{
  \item{pheno}{
Data frame where the first column is the line name (gid). The remaining columns can be either a phenotype or the levels of a fixed effect.  Any column not designated as a fixed effect is assumed to be a phenotype.
}
  \item{geno}{
Data frame with the marker names in the first column.  The second and third columns contain the chromosome and map position (either bp or cM), respectively, which are used only when plot=TRUE to make Manhattan plots.  If the markers are unmapped, just use a placeholder for those two columns.  Columns 4 and higher contain the marker scores for each line, coded as \{-1,0,1\} = \{aa,Aa,AA\}.  Fractional (imputed) and missing (NA) values are allowed.  The column names must match the line names in the "pheno" data frame.
}
\item{fixed}{
An array of strings containing the names of the columns that should be included as (categorical) fixed effects in the mixed model.
}
\item{K}{
Kinship matrix for the covariance between lines due to a polygenic effect.  If not passed, it is calculated from the markers using \code{\link{A.mat}}.
}
\item{n.PC}{
Number of principal components to include as fixed effects.  Default is 0 (equals K model).
}
\item{min.MAF}{
Specifies the minimum minor allele frequency (MAF).  If a marker has a MAF less than min.MAF, it is assigned a zero score.  
}
\item{n.core}{
Setting n.core > 1 will enable parallel execution on a machine with multiple cores (use only at UNIX command line).  
}
\item{P3D}{
When P3D=TRUE, variance components are estimated by REML only once, without any markers in the model.  When P3D=FALSE, variance components are estimated by REML for each marker separately. 
}
\item{plot}{
When plot=TRUE, qq and Manhattan plots are generated.
}
}
\details{
For unbalanced designs where phenotypes come from different environments, the environment mean can be modeled using the fixed option (e.g., fixed="env" if the column in the pheno data.frame is called "env").  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix.

The terminology "P3D" (population parameters previously determined) was introduced by Zhang et al. (2010).  When P3D=FALSE, this function is equivalent to EMMA with REML (Kang et al. 2008).  When P3D=TRUE, it is equivalent to EMMAX (Kang et al. 2010).  The P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

The dashed line in the Manhattan plots corresponds to an FDR rate of 0.05 and is calculated using the qvalue package (Storey and Tibshirani 2003).  The p-value corresponding to a q-value of 0.05 is determined by interpolation.  When there are no q-values less than 0.05, the dashed line is omitted.  
}
\value{
Returns a data frame where the first three columns are the marker name, chromosome, and position, and subsequent columns are the marker scores \eqn{(-log_{10}p)} for the traits.
}
\references{
Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Kang et al. 2010. Variance component model to account for sample structure in genome-wide association studies.
Nat. Genet. 42:348-354.

Storey and Tibshirani. 2003. Statistical significance for genome-wide studies. PNAS 100:9440-9445.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208. 

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide association studies. Nat. Genet. 42:355-360.
}

\examples{
#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
colnames(M) <- 1:200
geno <- data.frame(marker=1:1000,chrom=rep(1,1000),pos=1:1000,M,check.names=FALSE)

QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(M,u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

pheno <- data.frame(line=1:200,y=y)
scores <- GWAS(pheno,geno,plot=FALSE)
}
