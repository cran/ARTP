% --- Source file: ARTP.Rd ---
\name{ARTP}
\alias{ARTP}
\docType{package}
\title{
 Gene and pathway p-values using the Adaptive Rank Truncated Product test
}
\description{
An R package for computing gene and pathway p-values using the Adaptive Rank Truncated test.
This package can be used to analyze pathways/genes based on a genetic association study, with
either a continuous or a binary case-control outcome.
}
\details{ It is increasingly recognized that pathway analyses-a joint test of association between the outcome and a group of 
single nucleotide polymorphisms (SNPs) within a biological pathway-could potentially complement single-SNP analysis and
 provide additional insights for the genetic architecture of complex diseases. 
Building upon existing P-value combining methods, we propose a class of highly flexible pathway analysis approaches based 
on an adaptive rank truncated product statistic that can effectively combine evidence of associations over different SNPs and 
genes within a pathway. 
The statistical significance of the pathway-level test statistics is evaluated using a highly efficient permutation 
algorithm that remains computationally feasible irrespective of the size of the pathway and complexity of the underlying test
 statistics for summarizing SNP- and gene-level associations. 
We demonstrate through simulation studies that a gene-based analysis that treats the underlying genes, as opposed to the underlying SNPs,
 as the basic units for hypothesis testing, is a very robust and powerful approach to pathway-based association testing. \cr

 The function \code{\link{ARTP_pathway}} is used to compute gene and pathway p-values provided that the observed and 
 permutation p-values for each SNP already exist in files. The input files required
  for \code{\link{ARTP_pathway}} can be obtained by calling the function \code{\link{runPermutations}}.
}
\references{ 
Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, Caporaso N, Kraft P, Chatterjee N Pathway analysis by adaptive combination of P-values
 Genet Epidemiol 33(8):700-9; 2009 Dec.
}
\author{Kai Yu <yuka@mail.nih.gov> and William Wheeler <wheelerb@imsweb.com>}
\keyword{package}