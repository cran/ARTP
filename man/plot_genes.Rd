% --- Source file: plot_genes.Rd ---
\name{plot_genes}
\alias{plot_genes}
\title{Gene Plot}
\description{
  Plot the observed SNP p-values for each gene
}
\usage{
plot_genes(obs.outfile, gene.list, op=NULL)
}
\arguments{
  \item{obs.outfile}{The output file of observed p-values from \code{\link{runPermutations}}}
  \item{gene.list}{See \code{\link{gene.list}}}
  \item{op}{List of options. See \code{details}.}
}
\value{
  A data frame containing the SNP ids, parameter estimates, genes, etc.
 }
\details{ 
  If the option \code{gene.list$chrm.var} is not specified, then it is assumed that all the SNPs
  are on the same chromosome, and the same color will be used in the plot. If \code{gene.list$chrm.var}
  is specified, then the genes will be grouped by chromosome with the same color for the group.

  \bold{Options list:} \cr
 Below are the names for the options list \code{op}. 
\itemize{
  \item{\code{cex.axis}}{ See \code{\link[graphics]{par}} }
  \item{\code{colors}}{ Colors to use for each gene in the plot }
  \item{\code{maxLabelLen}}{ Maximum length of x-axis labels }
  \item{\code{chr.text}}{ See \code{\link[graphics]{par}} }
  \item{\code{x.las}}{ See \code{\link[graphics]{par}} }
  \item{\code{x.padj}}{ See \code{\link[graphics]{par}} }    

} % END: itemize



} % END: details

%\references{ }
%\author{Kai Yu }
\seealso{ \code{\link{runPermutations}} }
\examples{
 # Get the file of observed p-values
 obs_file <- system.file("sampleData", "obs_pvalues.txt", package="ARTP")

 # Define the gene-SNP list
 gs_file <- system.file("sampleData", "gene_SNP_data.txt", package="ARTP") 
 gene.list <- list(file=gs_file, delimiter="\t", header=1, 
                   snp.var="SNP", gene.var="Gene")
 
 plot_genes(obs_file, gene.list) 
  
}

\keyword{misc}

