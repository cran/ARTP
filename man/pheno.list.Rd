% --- Source file: pheno.list.Rd ---
\name{pheno.list}
\alias{pheno.list}
\docType{data}
\title{
 List to describe the covariate and outcome data 
}
\description{
 The list to describe the covariate and outcome data for \code{\link{runPermutations}}
}
\format{
The format is a list:
\describe{
   \item{file}{ Covariate data file. This file must have variable names, two
                 of which being an id variable and a response variable (see 
               \code{id.var} and \code{response.var}).
                  No default.}
   \item{id.var}{ Name of the id variable. 
                  No default.}
   \item{response.var}{ Name of the response variable. For logistic regression analyses, this variable
                 must be coded as 0 (control) and 1 (case).
                 No default. }
   \item{main.vars}{Character vector of variables names for variables in \code{file} 
                 that will be included in the model as main effects.
                 The default is NULL.}
   \item{delimiter}{ The delimiter in \code{file}.
                  The default is "".}
   \item{in.miss}{ Vector of character strings to define the missing values. This option
                   corresponds to the option \code{na.strings} in \code{\link[utils]{read.table}}.
                   The default is "NA".}
  }
}
\details{ In this list, \code{file}, \code{id.var}, and \code{response.var} must be specified.
   The variable \code{id.var} is the link between the covariate data and the genotype data.
   For each subject id, there must be the same subject id in the genotype data for that subject
   to be included in tha analysis. \cr
  \bold{Missing data:} If any of the variables defined in \code{main.vars}, \code{int.vars},
    or \code{response.var} contain missing values, then those subjects
  will be removed from the covariate and outcome data. After the subjects with missing values
  are removed, the subject ids are matched with the genotype data. 
}
%\examples{ Suppose the data of covariates is the space delimited file c:/temp/data/covData.txt 
%    with id variable "subjectID", and case-control variable "y".
% pheno.list <- list(file="c:/temp/data/covData.txt", file.type=3, id.var="subjectID",
%                response.var="y", delimiter=" ") 
%}
\keyword{misc}


