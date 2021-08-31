#'Summary method for dearseq objects 
#'
#'@param object an object of class \code{dear_seq}
#'
#'@param signif_threshold a value between \code{0} and {1} specifying the 
#'nominal significance threshold. Default is \code{0.05}.
#'
#'@param ... further arguments
#'
#'
#'@author Boris Hejblum
#'
#'@return a \code{list}
#'
#'@export

summary.dearseq <- function(object, signif_threshold = 0.05, ...){
  adj_pval <- object$pvals$adjPval
  names(adj_pval) <- rownames(object$pvals)
  
  adj_pval <- sort(adj_pval)
  
  s <- list()
  s$which_test <- object$which_test
  s$preprocessed <- object$preprocessed
  s$n_tests <- length(adj_pval)
  s$n_signif <- sum(adj_pval < signif_threshold)
  s$which_signif <- names(adj_pval)[adj_pval < signif_threshold]
  s$signif_threshold <- signif_threshold
  s$adj_pval <- adj_pval[adj_pval < signif_threshold]
  s$gsa <- !is.null(object$genesets)
  
  class(s) <- "summary.dearseq"
  return(s)
}

#'@rdname summary.dearseq
#'
#'@param x an object of class '\code{summary.dearseq}'.
#'
#'@export
print.summary.dearseq <- function(x, ...){
  if(x$gsa){
    cat(paste("dearseq identifies", x$n_signif, "significant genes sets out of", x$n_tests, "\n"))
  }else{
    cat(paste("dearseq identifies", x$n_signif, "significant genes out of", x$n_tests, "\n"))
  }
  cat(paste(" - test:", x$which_test, "\n"))
  cat(paste(" - FDR significance level:", x$signif_threshold, "\n"))
}
