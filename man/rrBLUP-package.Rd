\name{rrBLUP-package}
\alias{rrBLUP}
\docType{package}
\title{Ridge regression and other kernels for genomic selection}
\description{This package has been developed primarily for genomic prediction with mixed models (but it can also do genome-wide association mapping with \code{\link{GWAS}}). The heart of the package is the function \code{\link{mixed.solve}}, which is a general-purpose solver for mixed models with a single variance component other than the error.  Genomic predictions can be made by estimating marker effects (RR-BLUP) or by estimating line effects (G-BLUP). In Endelman (2011) I made the poor choice of using the letter G to denotype the genotype or marker data.  To be consistent with Endelman (2011) I have retained this notation in \code{\link{kinship.BLUP}}.  However, that function has now been superseded by \code{\link{kin.blup}} and \code{\link{A.mat}}, the latter being a utility for estimating the additive relationship matrix (A) from markers.  In these newer functions I adopt the usual convention that G is the genetic covariance (not the marker data), which is also consistent with the notation in Endelman and Jannink (2012).
}
\section{Parallel computing}{
For Mac, Linux, and UNIX users, R package multicore can be used in conjunction with rrBLUP to take advantage of multiple processors on a single machine.  This is useful for Gaussian kernel predictions and when using the EM imputation algorithm in \code{\link{A.mat}}.  You need R >= 2.14.1 for this to work properly, and you must also use R from the command line (not the GUI).  
}
\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. doi: 10.3835/plantgenome2011.08.0024

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics 2:1405-1413. doi: 10.1534/g3.112.004259
}
