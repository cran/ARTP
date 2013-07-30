% --- Source file: runPermutations.Rd ---
\name{runPermutations}
\alias{runPermutations}
\title{Calculate observed and permutation p-values for SNPs}
\description{
  Calculate observed and permutation p-values for SNPs
}
\usage{
runPermutations(snp.list, pheno.list, family, op=NULL)
}
\arguments{
  \item{snp.list}{A list describing the SNP data. See \code{\link{snp.list}} }
  \item{pheno.list}{A list describing the covariate and response data. See \code{\link{pheno.list}}}
  \item{family}{ 1 or 2,  1 = logistic regression,  2 = linear regression.}  
  \item{op}{List of options. See \code{details}.}
}
\value{
  The returned value is NULL, however 2 output files are created as defined by \code{op$obs.outfile}
  and \code{op$perm.outfile}.
 }
\details{  
  This function first reads the data stored in the files defined by \code{\link{snp.list}} and
  \code{\link{pheno.list}}.
  The subject ids in \code{snp.list$file} and \code{pheno.list$file} are matched and any subject not in both files will
  be removed from the analysis. Also, any subject with a missing value for the response or covariate
  will be removed.
  The function \code{\link{single.marker.test}} is called for each observed SNP and for each permutation.
  Depending on the response variable and genotype frequency counts, \code{\link{single.marker.test}} will call 
  \code{\link[stats]{glm.fit}}, \code{\link[stats]{fisher.test}} or \code{\link[stats]{lm}}. \cr
  Running a large number of permutations on a single machine could take a long time. However, if the user
  has access to multiple machines, then the permutations can be broken up across the different machines for
  faster computation. A different random seed should be set for each machine, and the output permutation files would need 
  to be combined into a single file before calling \code{\link{ARTP_pathway}}.

  \bold{Options list:} \cr
 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item{\code{nperm}}{ Number of permutations. The default is 100.}
  \item{\code{obs.outfile}}{ Output file for the observed results. The default is "obs.txt". }
  \item{\code{perm.outfile}}{ Output file for the permuted results. The default is "perm.txt".}
  \item{\code{perm.method}}{ 1 or 2 for the type of permutation. 1 is to permute the SNPs.
                            2 is to generate a new response using the base model.
                            For a continuous response, the residuals from the base model are permuted
                            and then added to the linear predictors from the base model to give 
                            the new response vector. For a binary response, the new response vector is
                            rbinom(n, 1, vals), where vals are the fitted values from the base model.
                            The default is 2. }
  \item{\code{min.count}}{ See \code{\link{single.marker.test}}. The default is 5. }
  \item{\code{miss.rate}}{ Maximum missing rate to include SNPs. Any SNP with missing rate greater than
                           \code{miss.rate} will be excluded. The default is 0.20. }  

} % END: itemize

\code{obs.outfile} will be a comma delimited file containing 5 rows: \cr
Row 1 contains the SNP ids. \cr
Row 2 contains the SNP p-values. \cr
Row 3 contains a value for how the p-value was computed (see the details of \code{\link{single.marker.test}}). \cr
Row 4 contains the estimate of the SNP main effect. \cr
Row 5 contains the estimated standard error of the SNP main effect. \cr
\code{perm.outfile} will be a comma delimited file, where each row are the permutation p-values
                   for all SNPs.

} % END: details

%\references{ }
\author{Kai Yu and William Wheeler}
\seealso{ \code{\link{single.marker.test}} \code{\link{snp.list}} \code{\link{pheno.list}} }
\examples{
  # Define snp.list
 geno_file <- system.file("sampleData", "geno_data.txt", package="ARTP")
 snp.list <- list(file=geno_file, file.type=2, delimiter="\t")

 # Define pheno.list
 pheno_file <- system.file("sampleData", "pheno_data.txt", package="ARTP")
 pheno.list <- list(file=pheno_file, delimiter="\t", id.var="ID",
                    response.var="Y", main.vars=c("X1", "X2"))

 # Options list. Change obs.outfile and perm.outfile if needed.
 op <- list(nperm=10, obs.outfile="./obs.txt", perm.outfile="./perm.txt",
            perm.method=2)

 # Not run
 # runPermutations(snp.list, pheno.list, 1, op=op)  

}

\keyword{model} 

