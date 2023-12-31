% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vc_test_perm.R
\name{vc_test_perm}
\alias{vc_test_perm}
\title{Permutation-based variance component test statistic and p-value}
\usage{
vc_test_perm(
  y,
  x,
  indiv = rep(1, nrow(x)),
  phi,
  w,
  Sigma_xi = diag(ncol(phi)),
  n_perm = 1000,
  progressbar = TRUE,
  parallel_comp = TRUE,
  nb_cores = parallel::detectCores(logical = FALSE) - 1,
  genewise_pvals = FALSE,
  adaptive = TRUE,
  max_adaptive = 64000,
  homogen_traj = FALSE,
  na.rm = FALSE
)
}
\arguments{
\item{y}{a numeric matrix of dim \code{G x n} containing the raw RNA-seq
counts for G genes from \code{n} samples.}

\item{x}{a numeric design matrix of dim \code{n x p} containing the \code{p}
covariates to be adjusted for.}

\item{indiv}{a vector of length \code{n} containing the information for
attributing each sample to one of the studied individuals. Coerced to be a
\code{factor}.}

\item{phi}{a numeric design matrix of size \code{n x K} containing the
\code{K} variables to be tested}

\item{w}{a vector of length \code{n} containing the weights for the \code{n}
samples.}

\item{Sigma_xi}{a matrix of size \code{K x K} containing the covariance matrix
of the \code{K} random effects.}

\item{n_perm}{the number of perturbations. Default is \code{1000}.}

\item{progressbar}{logical indicating wether a progressBar should be displayed
when computing permutations (only in interactive mode).}

\item{parallel_comp}{a logical flag indicating whether parallel computation
should be enabled. Only Linux and MacOS are supported, this is ignored on
Windows. Default is \code{TRUE}.}

\item{nb_cores}{an integer indicating the number of cores to be used when
\code{parallel_comp} is \code{TRUE}.
Default is \code{parallel::detectCores(logical=FALSE) - 1}.}

\item{genewise_pvals}{a logical flag indicating whether gene-wise p-values
should be returned. Default is \code{FALSE} in which case gene-set p-value is
computed and returned instead.}

\item{adaptive}{a logical flag indicating whether adaptive permutation should
be performed. Default is \code{TRUE}}

\item{max_adaptive}{The maximum number of permutations considered.
Default is \code{64000}}

\item{homogen_traj}{a logical flag indicating whether trajectories should be
considered homogeneous. Default is \code{FALSE} in which case trajectories are
not only tested for trend, but also for heterogeneity.}

\item{na.rm}{logical: should missing values (including \code{NA} and
\code{NaN}) be omitted from the calculations? Default is \code{FALSE}.}
}
\value{
A list with the following elements when the set p-value is computed:
\itemize{
  \item \code{set_score_obs}: the approximation of the observed set score
  \item \code{set_pval}: the associated set p-value
}
or a list with the following elements when gene-wise p-values are computed:
\itemize{
  \item \code{gene_scores_obs}: vector of approximating the observed
  gene-wise scores
  \item \code{gene_pvals}: vector of associated gene-wise p-values
  \item \code{ds_fdr}: vector of associated gene-wise discrete false
  discovery rates
}
}
\description{
This function computes an approximation of the Variance Component test for a
mixture of \eqn{\chi^{2}}s using permutations. This is preferable to the
asymptotic approximation for small sample sizes. We rely on exact p-values
following Phipson and Smyth, 2010 (see References).
}
\examples{
set.seed(123)

##generate some fake data
########################
n <- 100
r <- 12
t <- matrix(rep(1:3), 4, ncol=1, nrow=r)
sigma <- 0.4
b0 <- 1

#under the null:
b1 <- 0
#under the alternative:
b1 <- 0.5
y.tilde <- b0 + b1*t + rnorm(r, sd = sigma)
y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
      matrix(rep(y.tilde, n), ncol=n, nrow=r))
x <- matrix(1, ncol=1, nrow=r)

#run test
permTestRes <- vc_test_perm(y, x, phi=t,
                           w=matrix(1, ncol=ncol(y), nrow=nrow(y)),
                           indiv=rep(1:4, each=3), n_perm=50, #1000,
                           parallel_comp = FALSE)
permTestRes$set_pval

}
\references{
Phipson B, and Smyth GK (2010). Permutation p-values should
never be zero: calculating exact p-values when permutations are randomly
drawn. \emph{Statistical Applications in Genetics and Molecular Biology},
Volume 9, Issue 1, Article 39.
\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}
}
