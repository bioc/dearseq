#'Plotting raw p-values histogram
#'
#'Display the histogram of raw p-values for diagnostic plots
#'
#'@param pvals a vector of raw p-values
#'
#'@param binwidth a value specifying the width of the histogram bins. 
#'Default is \code{0.02}.
#'
#'@author Boris Hejblum
#'
#'@return a \code{\link[ggplot2]{ggplot}} object
#'
#'@import ggplot2
#'@export
#'
#'@examples
#'#generate fake data
#'n <- 1000 #number of genes
#'nr=5 #number of measurements per subject (grouped data)
#'ni=50 #number of subjects
#'r <- nr*ni #number of measurements
#'t <- matrix(rep(1:nr), ni, ncol=1, nrow=r) # the variable to be tested 
#'sigma <- 0.5
#'y.tilde <- rnorm(r, sd = sigma)
#'y <- t(matrix(rnorm(n*r, sd = sqrt(sigma*abs(y.tilde))), ncol=n, nrow=r) +
#'       matrix(rep(y.tilde, n), ncol=n, nrow=r))
#'
#'#run test
#'#asymptotic test with preprocessed grouped data
#'res_genes <- dear_seq(exprmat=y, covariates=x, variables2test=t,
#'                    sample_group=rep(1:ni, each=nr),
#'                    which_test = "asymptotic",
#'                    which_weights='none', preprocessed=TRUE)
#'plot_hist_pvals(res_genes$pvals$rawPval)
#'
plot_hist_pvals <- function(pvals, binwidth = 0.02){
  
  stopifnot(is.numeric(pvals))
  
  df2plot <- cbind.data.frame("Rawpvals" = pvals)
  
  ggp <- ggplot(df2plot) + 
    geom_histogram(aes_string(x = "Rawpvals"), color="white", binwidth = binwidth) +
    xlab("Raw p-values") +
    theme_bw() +
    ggtitle("Histogram of raw p-values", subtitle = "before multiple testing correction")
  
  return(ggp)
}
