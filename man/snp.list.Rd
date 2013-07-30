% --- Source file: snp.list.Rd ---
\name{snp.list}
\alias{snp.list}
\docType{data}
\title{
 List to describe the genotype data
}
\description{
 The list to describe the genotype data for \code{\link{runPermutations}}
}
\format{
  The format is a list:
  \describe{
  \item{file}{File to use. No default. }
  \item{file.type}{ 2 or 3 (see details). }
  \item{delimiter}{ The delimiter used in \code{file}. }
  \item{in.miss}{Vector of values to denote the missing values in \code{file}.
                 The default is "  " (2 blank spaces). }
  \item{heter.codes}{ Vector of codes used for the heterozygous genotype.
                    If NULL, then it is assumed that the heterozygous genotype 
                    is of the form "AB", "Aa", "CT", ... etc, ie a 2-character
                    string with different characters (case sensitive).
                    The default is NULL.}
  \item{id.var}{   (Only for \code{file.type} = 3) The subject id variable.
                    The default is 1.}
  }
}

\details{ In this list, \code{file} must be specified. 
         If the SNPs are coded in the standard (0,1,2) coding, then
         set \code{heter.codes} to 1 (the heterozygous genotype).
          \cr \cr
                     \code{Type 2} has data in the form: \cr \cr
                   \tabular{cccc}{
                       \tab subject1 \tab subject2 \tab subject3 \cr
                      snp1 \tab 0 \tab 2 \tab 1 \cr
                      snp2 \tab 1 \tab 1 \tab 0 
                    }
                    The first row
                    must contain the subject ids. Starting from row 2, the first
                    delimited field must contain the SNP id. The remaining
                    delimited fields contain the genotypes.
                    Rows are SNPs, columns are the subjects. \cr
      
                    \code{Type 3} has data of the form: \cr \cr
                    \tabular{ccc}{
                      id \tab snp1 \tab snp2  \cr
                      subject1 \tab 0 \tab 1 \cr
                      subject2 \tab 2 \tab 1 \cr
                      subject3 \tab 1 \tab 0 
                    }
                   
}
\examples{
# Suppose the genotype data is a tab-delimited, type 2 file: c:/temp/data/geno1.txt.
#  Also assume the data has the trend coding 0, 1, 2 with NA as missing values.
# The below list is for processing the file.
snp.list <- list(file="C:/temp/data/geno1.txt", delimiter="\t", file.type=2, 
                 heter.codes=1, in.miss=NA)
}
\keyword{misc}

