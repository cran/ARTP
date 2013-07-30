% --- Source file: gene.list.Rd ---
\name{gene.list}
\alias{gene.list}
%\docType{package}
\title{
 List to describe the gene-SNP file
}
\description{
 The list to describe the gene-SNP file for \code{\link{ARTP_pathway}}
}
\format{
 The format is a list:
\describe{
  \item{file}{Text file containing at least 2 columns, where one column is for the SNPs and
              the other column is for the gene containing the SNP. No default. }
  \item{delimiter}{ The delimiter used in \code{file}. }
  \item{gene.var}{ Variable name or column number of the gene variable. The default is "Gene".}
  \item{snp.var}{ Variable name or column number of the SNP variable. The default is "SNP".}
  \item{chrm.var}{ Variable name or column number of the chromosome variable. The default is "Chr".
                  This option is only used for calling \code{\link{plot_genes}}.}
  \item{header}{ 0 or 1 to denote if \code{file} contains a header of variable names. The default is 1.}
  }
}
\details{ 
  All the genes and SNPs listed in this file define a single pathway.
}
%\examples{}
\keyword{misc}

