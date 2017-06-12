\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates the realized additive relationship matrix.
}
\usage{
A.mat(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
      n.core=1,shrink=FALSE,return.imputed=FALSE)
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency. The A matrix is not sensitive to rare alleles, so by default only monomorphic markers are removed.
}
\item{max.missing}{
Maximum proportion of missing data; default removes completely missing markers.
}
\item{impute.method}{
There are two options. The default is "mean", which imputes with the mean for each marker.  The "EM" option imputes with an EM algorithm (see details).
}
\item{tol}{
Specifies the convergence criterion for the EM algorithm (see details). 
}
\item{n.core}{
Specifies the number of cores to use for parallel execution of the EM algorithm (use only at UNIX command line).
}
\item{shrink}{
Set shrink=FALSE to disable shrinkage estimation.  See Details for how to enable shrinkage estimation. 
}
\item{return.imputed}{
When TRUE, the imputed marker matrix is returned.  
}
}
\details{
At high marker density, the relationship matrix is estimated as \eqn{A=W W'/c}, where \eqn{W_{ik} = X_{ik} + 1 - 2 p_k} and \eqn{p_k} is the frequency of the 1 allele at marker k.  By using a normalization constant of \eqn{c = 2 \sum_k {p_k (1-p_k)}}, the mean of the diagonal elements is \eqn{1 + f} (Endelman and Jannink 2012).

The EM imputation algorithm is based on the multivariate normal distribution and was designed for use with GBS (genotyping-by-sequencing) markers, which tend to be high density but with lots of missing data.  Details are given in Poland et al. (2012). The EM algorithm stops at iteration \eqn{t} when the RMS error = \eqn{n^{-1} \|A_{t} - A_{t-1}\|_2} < tol.  

Shrinkage estimation can improve the accuracy of genome-wide marker-assisted selection, particularly at low marker density (Endelman and Jannink 2012). The shrinkage intensity ranges from 0 (no shrinkage) to 1 (\eqn{A=(1+f)I}). Two algorithms for estimating the shrinkage intensity are available.  The first is the method described in Endelman and Jannink (2012) and is specified by \code{shrink=list(method="EJ")}.  The second involves designating a random sample of the markers as simulated QTL and then regressing the A matrix based on the QTL against the A matrix based on the remaining markers (Yang et al. 2010; Mueller et al. 2015). The regression method is specified by \code{shrink=list(method="REG",n.qtl=100,n.rep=5)}, where the parameters \code{n.qtl} and \code{n.rep} can be varied to adjust the number of simulated QTL and number of replicates, respectively.  

The shrinkage and EM-imputation options are designed for opposite scenarios (low vs. high density) and cannot be used simultaneously. 
When the EM algorithm is used, the imputed alleles can lie outside the interval [-1,1].  Polymorphic markers that do not meet the min.MAF and max.missing criteria are not imputed.
}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$A}{the A matrix}
\item{$imputed}{the imputed marker matrix}
}
}

\references{
Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics. 2:1405-1413. doi: 10.1534/g3.112.004259

Mueller et al. 2015. Shrinkage estimation of the genomic relationship matrix can improve genomic estimated breeding values in the training set. Theor Appl Genet doi: 10.1007/s00122-015-2464-6

Poland, J., J. Endelman et al. 2012. Genomic selection in wheat breeding using genotyping-by-sequencing. Plant Genome 5:103-113. doi: 10.3835/plantgenome2012.06.0006

Yang et al. 2010. Common SNPs explain a large proportion of the heritability for human height. Nat. Genetics 42:565-569.
}

\examples{
#random population of 200 lines with 1000 markers
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

A <- A.mat(X)
}