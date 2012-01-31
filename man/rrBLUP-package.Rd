\name{rrBLUP-package}
\alias{rrBLUP}
\docType{package}
\title{Ridge regression and other kernels for genomic selection}
\description{This package has been designed for both genomic selection and association mapping.  Some details of using the package for genomic selection have been published in the Plant Genome.  The heart of the package is the function \code{\link{mixed.solve}}, which can be used to model marker effects as random effects or the genotypic values of the lines as random effects.  In the latter case, the function \code{\link{A.mat}} is useful for calculating the additive relationship matrix and the prediction of breeding values.  To include epistatic effects in the genotypic value predictions, use the function \code{\link{kinship.BLUP}}.   
}
\section{Association mapping}{
Use function \code{\link{GWA}} for association mapping.
}
\section{Missing data}{
A number of improvements have been made concerning the handling of missing data since the original publication of the package.  The functions \code{\link{mixed.solve}}, \code{\link{kinship.BLUP}}, and \code{\link{GWA}} will tolerate missing phenotypic data: those observations are simply not used.  When genotypic data are missing, both \code{\link{kinship.BLUP}} (option "RR") and \code{\link{GWA}} rely on the EM algorithm in \code{\link{A.mat}}, which can also be used in conjunction with \code{\link{mixed.solve}}. The non-additive kernels in \code{\link{kinship.BLUP}} are based on \code{\link{dist}}, which uses pairwise complete observations.  
}
\section{Parallel computing}{
For Mac, Linux, and UNIX users, R package multicore can be used in conjunction with rrBLUP to take advantage of multiple processors on a single machine.  For large data sets, especially when there is missing data, I highly recommend trying out this feature, which is available with \code{\link{kinship.BLUP}}, \code{\link{A.mat}}, and \code{\link{GWA}}.  You need R >= 2.14.1 for this to work properly, and you must also use R from the command line (not the GUI).  
}
\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255.
}
