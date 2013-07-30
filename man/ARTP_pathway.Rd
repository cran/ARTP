% --- Source file: ARTP_pathway.Rd ---
\name{ARTP_pathway}
\alias{ARTP_pathway}
\title{Gene and pathway p-values using ARTP}
\description{
  Calculate gene and pathway p-values using the Adaptive Rank Truncated Product test
}
\usage{
ARTP_pathway(obs.file, perm.file, nperm, temp.dir, gene.list=NULL, op=NULL)
}
\arguments{
  \item{obs.file}{The output file \code{obs.outfile} from \code{\link{runPermutations}} or
                  a file with the SNP ids and p-values (see details). }
  \item{perm.file}{The output file \code{perm.outfile} from \code{\link{runPermutations}} or a files with the SNP ids
                   and p-values (see details).}
  \item{nperm}{The number of permutations in the output file \code{perm.outfile} from \code{\link{runPermutations}} }
  \item{temp.dir}{A folder to keep temporary files that will be created. }
  \item{gene.list}{A list describing the gene-SNP data. See \code{\link{gene.list}}. 
                   If NULL, then it is assumed that all SNPs belong to the same gene.
                   The default value is NULL.}
  \item{op}{List of options. See \code{details}.}
}
\value{
  The returned value is a list with names "pathway.pvalue" and "gene.table".
  pathway.pvalue is the ARTP p-value for the pathway. gene.table is a data frame containing the 
  gene name, number of SNPs in the gene that were included in the analysis, and the ARTP p-value
  for the gene.
 }
\details{  
   If the p-values are not computed using \code{\link{runPermutations}}, then the format for \code{obs.file} and
    \code{perm.file} should be as follows. Both files must be uncompressed, comma seperated files with the first row as
   the SNP ids in the same order. Row 2 of \code{obs.file} has the observed p-values, and starting from row 2 in 
   \code{perm.file} are the permuted p-values. 

   A random seed should be set before calling \code{ARTP_pathway} in order to reproduce results. 
   The randomness is due to the ranking of p-values, where ties are broken randomly. \cr 

  \bold{Options list:} \cr
 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item{\code{inspect.snp.n}}{ The number of candidate truncation points to inspect the top SNPs in a gene.
                               The default is 1. }
  \item{\code{inspect.snp.percent}}{ A value \code{x} between 0 and 1 such that a truncation point will be defined
                                    at every \code{x} percent of the top SNPs. 
                                    The default is 0 so that the truncation points will be 1:\code{inspect.snp.n}.}
  \item{\code{inspect.gene.n}}{ The number of candidate truncation points to inspect the top genes
                                in the pathway. The default is 10. }
  \item{\code{inspect.gene.percent}}{ A value \code{x} between 0 and 1 such that a truncation point will be defined
                                    at every \code{x} percent of the top genes. The default is 0.05. }
  %\item{\code{pathway.method}}{ 0 or 1 for the method to compute the pathway p-value.
  %                              Method 0 uses the Adaptive Rank Truncated Product test. 
  %                              Method 1 uses prior probabilities for each gene.
  %                              The default is 0. }
  %\item{\code{parms.list}}{ List to specify the parameters for \code{pathway.method} = 1. This list
  %                         can have the names "theta", "file", "delimiter", "header", "gene.var", "parm.var".
  %                         The default value of \code{theta} is 0.2. If \code{file} is not specified,
  %                         then the \code{pi} parameter for each gene is set to 1. }

} % END: itemize

Assume the number of SNPs in a gene is 100. Below are examples of the truncation points for different
 values of \code{inspect.snp.n} and \code{inspect.snp.percent}. 

\tabular{rrr}{
inspect.snp.n \tab inspect.snp.percent \tab truncation points \cr
1 \tab 0 \tab 1 \cr
1 \tab 0.05 \tab 5 \cr
1 \tab 0.25 \tab 25 \cr
1 \tab 1 \tab 100 \cr
2 \tab 0 \tab 1, 2 \cr
2 \tab 0.05 \tab 5, 10 \cr
2 \tab 0.25 \tab 25, 50 \cr
2 \tab 1 \tab 100 \cr
3 \tab 0.2 \tab 20, 40, 60 
} % END: tabular


} % END: details

%\references{ }
\author{Kai Yu }
\seealso{ \code{\link{runPermutations}} }
\examples{
  # Get the file of observed p-values
 obs_file <- system.file("sampleData", "obs_pvalues.txt", package="ARTP")

 # Get the file of permutation p-values
 perm_file <- system.file("sampleData", "perm_pvalues.txt", package="ARTP")

 # Define the gene-SNP list
 gs_file <- system.file("sampleData", "gene_SNP_data.txt", package="ARTP") 
 gene.list <- list(file=gs_file, delimiter="\t", header=1, 
                   snp.var="SNP", gene.var="Gene")

 # Call the ARTP function
 nperm <- 100  # The number of permutations in perm_file
 temp.dir <- "C:/temp/"
 set.seed(123)
 # ARTP_pathway(obs_file, perm_file, nperm, temp.dir, gene.list=gene.list)

 # Now asume that all SNPs belong to the same gene
 # ARTP_pathway(obs_file, perm_file, nperm, temp.dir)

}

\keyword{model} 

