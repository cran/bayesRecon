% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrink_cov.R
\name{schaferStrimmer_cov}
\alias{schaferStrimmer_cov}
\title{Schäfer Strimmer covariance shrinkage}
\usage{
schaferStrimmer_cov(x)
}
\arguments{
\item{x}{matrix of samples with dimensions nxp (n samples, p dimensions).}
}
\value{
A list containing the shrinkage estimator and the optimal lambda. The list has the following named elements:
\itemize{
\item \code{shrink_cov}: the shrinked covariance matrix (\code{p} x \code{p});
\item \code{lambda_star}: the optimal lambda for the shrinkage;
}
}
\description{
Computes the Schäfer Strimmer shrinkage estimator for a covariance matrix
from a matrix of samples.
}
\details{
This function computes the shrinkage to a diagonal covariance with unequal variances.
Note that here we use the estimators \eqn{S = X X^T/n} and \eqn{T = diag(S)} and we internally
use the correlation matrix in place of the covariance to compute the optimal shrinkage factor.
}
\examples{

# Generate some multivariate normal samples
# Parameters
nSamples <- 200
pTrue <- 2

# True moments
trueSigma <- matrix(c(3,2,2,2), nrow=2)
chol_trueSigma <- chol(trueSigma)
trueMean <- c(0,0) 

# Generate samples
set.seed(42)
x <- replicate(nSamples, trueMean) +  
     t(chol_trueSigma)\%*\%matrix(stats::rnorm(pTrue*nSamples), 
                                nrow = pTrue, ncol = nSamples)
x <- t(x) 
res_shrinkage <- schaferStrimmer_cov(x)
res_shrinkage$lambda_star # should be 0.01287923

}
\references{
Schäfer, Juliane, and Korbinian Strimmer. (2005). \emph{A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics.} Statistical Applications in Genetics and Molecular Biology 4: Article32. \doi{10.2202/1544-6115.1175}.
}
