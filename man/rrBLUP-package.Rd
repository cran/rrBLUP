\name{rrBLUP-package}
\alias{rrBLUP}
\docType{package}
\title{Ridge regression and other kernels for genomic selection}
\description{This package has been developed primarily for genomic prediction with mixed models (but it can also do genome-wide association mapping with \code{\link{GWAS}}). The heart of the package is the function \code{\link{mixed.solve}}, which is a general-purpose solver for mixed models with a single variance component other than the error.  Genomic predictions can be made by estimating marker effects (RR-BLUP) or by estimating line effects (G-BLUP). In Endelman (2011) I made the poor choice of using the letter G to denotype the genotype or marker data.  To be consistent with Endelman (2011) I have retained this notation in \code{\link{kinship.BLUP}}.  However, that function has now been superseded by \code{\link{kin.blup}} and \code{\link{A.mat}}, the latter being a utility for estimating the additive relationship matrix (A) from markers.  In these newer functions I adopt the usual convention that G is the genetic covariance (not the marker data), which is also consistent with the notation in Endelman and Jannink (2012).
	
Vignettes illustrating some of the features of this package can be found at \url{http://potatobreeding.cals.wisc.edu/software}.
}

\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. doi: 10.3835/plantgenome2011.08.0024

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics 2:1405-1413. doi: 10.1534/g3.112.004259
}
