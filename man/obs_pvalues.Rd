% --- Source file: obs_pvalues.Rd ---
\name{obs_pvalues}
\alias{obs_pvalues}
%\docType{package}
\title{
 Observed p-values 
}
\description{
 Sample file of observed p-values for the example in \code{\link{ARTP_pathway}}
}
\details{ This is a comma delimited file where the first row contains the SNP ids, second row
 contains the p-values, third row contains the method of p-value computation (see \code{\link{single.marker.test}}),
 fourth row contains the estimated main efect of the SNP, and fifth row contains the estimated SNP main effect 
 standard error.
}
\examples{
  # Read in and print the data
  f <- system.file("sampleData", "obs_pvalues.txt", package="ARTP")
  x <- scan(f, what="character", sep="\n")
  substring(x, 1, 50)

}
\keyword{data}

