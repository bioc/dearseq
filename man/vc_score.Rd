% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vc_score.R
\name{vc_score}
\alias{vc_score}
\title{Computes variance component score test statistics}
\usage{
vc_score(y, x, indiv, phi, w, Sigma_xi = diag(ncol(phi)), na_rm = FALSE)
}
\arguments{
\item{y}{a numeric matrix of dim \code{g x n} containing the raw RNA-seq
counts for g genes from \code{n} samples.}

\item{x}{a numeric design matrix of dim \code{n x p} containing the \code{p}
covariates to be adjusted for.}

\item{indiv}{a vector of length \code{n} containing the information for
attributing each sample to one of the studied individuals. Coerced
to be a \code{factor}.}

\item{phi}{a numeric design matrix of size \code{n x K} containing the
\code{K} variables to .be tested}

\item{w}{a vector of length \code{n} containing the weights for the \code{n}
samples.}

\item{Sigma_xi}{a matrix of size \code{K x K} containing the covariance matrix
of the \code{K} random effects on \code{phi}.}

\item{na_rm}{logical: should missing values (including \code{NA} and
\code{NaN}) be omitted from the calculations? Default is \code{FALSE}.}
}
\value{
A list with the following elements:\itemize{
  \item \code{score}: approximation of the set observed score
  \item \code{q}: observation-level contributions to the score
  \item \code{q_ext}: pseudo-observations used to compute the covariance,
   taking into account the contributions of OLS estimates
  \item \code{gene_scores_unscaled}: a vector of the approximations of the
  individual gene scores
}
}
\description{
This function computes the variance component score test statistics
}
\examples{
set.seed(123)

##generate some fake data
########################
n <- 100
r <- 12
t <- matrix(rep(1:3), r/3, ncol=1, nrow=r)
sigma <- 0.4
b0 <- 1

#under the null:
b1 <- 0
#under the alternative:
b1 <- 0.7
y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
      matrix(rep(y.tilde, n), ncol=n, nrow=r))
x <- matrix(1, ncol=1, nrow=r)

#run test
scoreTest <- vc_score(y, x, phi=t, w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
                    Sigma_xi=matrix(1), indiv=rep(1:(r/3), each=3))
scoreTest$score

}
\keyword{internal}
