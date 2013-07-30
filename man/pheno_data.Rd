% --- Source file: pheno_data.Rd ---
\name{pheno_data}
\alias{pheno_data}
%\docType{package}
\title{
 Sample covariate and response data 
}
\description{
 Sample covariate and response data for \code{\link{runPermutations}}
}
\details{ The file pheno_data.txt is a
   tab-delimited type 3 data set (see \code{file.type} in \code{\link{pheno.list}}).
   It contains the variables:
 \itemize{
  \item{\code{ID}}{ The subject id}
  \item{\code{Y}}{ Case-control status (0, 1)}
  \item{\code{X1}}{ Continuous covariate}
  \item{\code{X2}}{ Continuous covariate}
 } 
}
\examples{
 # Load and print the first 5 rows 
 data(pheno_data, package="ARTP") 

 pheno_data[1:5, ]
}
\keyword{data}




