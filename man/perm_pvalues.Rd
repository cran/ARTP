% --- Source file: perm_pvalues.Rd ---
\name{perm_pvalues}
\alias{perm_pvalues}
%\docType{package}
\title{
 Permutation p-values 
}
\description{
 Sample file of permutation p-values for the example in \code{\link{ARTP_pathway}}
}
\details{ 
This is a comma delimited file with 101 rows with row 1 containing the SNP ids, and rows 2-101 
containing the p-values. Each row represents one permutation.
}
\examples{
# Read in and print the data
  f <- system.file("sampleData", "perm_pvalues.txt", package="ARTP")
  x <- scan(f, what="character", sep="\n")
  substring(x[1:5], 1, 50)
}
\keyword{data}

