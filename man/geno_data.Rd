% --- Source file: geno_data.Rd ---
\name{geno_data}
\alias{geno_data}
%\docType{package}
\title{
 Sample genotype data 
}
\description{
 Sample genotype data for \code{\link{runPermutations}}
}
\details{ geno_data.rda is a type 1 data file (see \code{file.type} in \code{\link{snp.list}}).
This data contains 50 SNPs and 500 subjects, and is tab delimited.
The first row of the data contains the subject ids. Starting from row 2, are the SNP ids
 and the genotypes for each subject. The genotypes are coded as AA, AG, GG.
}
\examples{
 # Load and print a substring the first 5 lines
 data(geno_data, package="ARTP")
 
 substring(geno_data[1:5], 1, 50)
}
\keyword{data}

