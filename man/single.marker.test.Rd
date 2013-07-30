% --- Source file: single.marker.test.Rd ---
\name{single.marker.test}
\alias{single.marker.test}
\title{Single SNP test}
\description{
  Perform an association test for 1 SNP
}
\usage{
single.marker.test(y, covariates, weights, offset, control, snpcol,
                   min.count=5, y.continuous=FALSE)
}
\arguments{
  \item{y}{The response vector}
  \item{covariates}{A design matrix where the SNP is the last column. The SNP must be coded as 0-1-2.}
  \item{weights}{Vector of weights.}
  \item{offset}{Vector for the offset.}
  \item{control}{List for \code{\link[stats]{glm.control}}}
  \item{snpcol}{Number of columns of the design matrix \code{covariates}}
  \item{min.count}{The minimum number of subjects to have in at least 2 of the genotype categories (0-1-2),
                   if \code{y} is continuous. If \code{y} is binary, then
                   the minimum expected frequency count for cases or controls to use logistic regression; 
                   otherwise, Fisher's exact test will be used. 
                   The default is 5}
  \item{y.continuous}{TRUE or FALSE for whether or not \code{y} is continuous. 
                      If FALSE, then \code{y} must be coded as 0-1.    
                        The default is FALSE.}
}
\value{
  The returned object is a vector of length 4 containing the p-value, p-value flag (see details),
      SNP main effect estimate, and standard error of the SNP main effect estimate. 
      If Fisher's exact test was used, then the main effect and standard error will be set to NA.
 }
\details{ 
 The input vectors and matrices must not contain missing values. 
 To compute the p-value, either \code{\link[stats]{glm.fit}}, \code{\link[stats]{fisher.test}} or
  \code{\link[stats]{lm}} is called.
 The p-value flag is a value for how the p-value was computed:
  \tabular{rrr}{
  Value \tab Genetic Model \tab Test \cr 
  0 \tab trend \tab Wald test from logistic/linear regression \cr
 -1 \tab dominant \tab Fisher's exact test \cr 
 -2 \tab recessive \tab Fisher's exact test \cr 
  1 \tab dominant \tab Wald test from logistic regression \cr
  2 \tab recessive \tab Wald test from logistic regression
  } % END: tabular


} % END: details

%\references{ }
\author{Kai Yu and Qizhai Li}
\seealso{ \code{\link{runPermutations}} }
\examples{

  # Generate data
  set.seed(123)
  n <- 1000
  y <- rbinom(n, 1, 0.5)
  snp <- rbinom(n, 2, 0.4)
  weights <- rep.int(1, times=n)
  offset  <- rep.int(0, times=n)
  control  <- glm.control()
 
  # Create a design matrix
  x <- matrix(data=NA, nrow=n, ncol=3)
  x[, 1] <- 1  # Intercept column
  x[, 2] <- runif(n) # Continuous covariate
  x[, 3] <- snp

  single.marker.test(y, x, weights, offset, control, 3)

}

\keyword{model} 

