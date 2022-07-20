#'Plotting mean-variance fit for precision weights estimation
#'
#'Display the variability with respect to the level of expression and the 
#'associated smoothed estimation of precision weights accounting for 
#'heteroscedasticity.
#'
#'@param x a list (such as outputed by the functions 
#'\code{\link[dearseq]{sp_weights}} or \code{\link[dearseq]{voom_weights}}) 
#'containing the following components:\itemize{
#'\item \code{weights}: a matrix \code{n x G} containing the estimated precision 
#'weights 
#'\item \code{plot_utilities}: a list containing the following elements:\itemize{
#'      \item\code{reverse_trans}: a function encoding the reverse function used 
#'      for smoothing the observations before computing the weights
#'      \item\code{method}: the weight computation method (either \code{"voom"} 
#'      or \code{"loclin"})
#'      \item\code{smth}: the vector of the smoothed values computed 
#'      \item\code{gene_based}: a logical indicating whether the computed 
#'      weights are based on average at the gene level or on individual 
#'      observations
#'      \item\code{mu}: the transformed observed counts or averages
#'      \item\code{v}: the observed variability estimates
#'}
#'}
#'
#'@author Boris Hejblum
#'
#'@return a \code{\link[ggplot2]{ggplot}} object
#'
#'@import ggplot2
#'@importFrom scattermore geom_scattermore
#'@export
#'
#'@examples
#'G <- 10000
#'n <- 12
#'p <- 2
#'y <- sapply(1:n, FUN = function(x){rnbinom(n = G, size = 0.07, mu = 200)})
#'x <- sapply(1:p, FUN = function(x){rnorm(n = n, mean = n, sd = 1)})
#'
#'if(interactive()){
#'  w <- sp_weights(y, x, use_phi=FALSE, na.rm = TRUE, gene_based=TRUE)
#'  plot_weights(w)
#'
#'  vw <-  voom_weights(y, x)
#'  plot_weights(vw)
#'}
#'
plot_weights <- function(x){
  
  stopifnot(!is.null(x$plot_utilities))
  stopifnot(!is.null(x$plot_utilities$method))
  
  if(x$plot_utilities$method == "voom"){
    r_tilde <- x$plot_utilities$mu
    s_rs <- x$plot_utilities$v
    o <- order(r_tilde, na.last = NA)
    
    plot_df <- data.frame(r_tilde_o = r_tilde[o], s_rs_o = s_rs[o])
    plot_df_lo <- data.frame(lo.x = x$plot_utilities$smth$x, lo.y = x$plot_utilities$smth$y)
    ggp <- (ggplot(data = plot_df) +
              geom_scattermore(aes_string(x = "r_tilde_o", y = "s_rs_o"),
                               alpha = 0.5, color = "grey30") +
              theme_bw() +
              xlab("Gene average") +
              ylab("sqrt(sd)") +
              ggtitle("Mean-variance") +
              geom_line(data = plot_df_lo,
                        aes_string(x = "lo.x", y = "lo.y"),
                        color = "blue", lwd = 1.4, lty = "solid",
                        alpha = 0.8))
    
    
  }else if(x$plot_utilities$method == "loclin"){
    
    mu_orig <- x$plot_utilities$reverse_trans(x$plot_utilities$mu)
    o <- order(mu_orig, na.last = NA)
    plot_df <- data.frame("m_o" = mu_orig[o], 
                          "v_o" = x$plot_utilities$v[o])
    smth <- x$plot_utilities$smth[o]
    if (x$plot_utilities$gene_based) {
      plot_df_lo <- data.frame("lo.x" = mu_orig[o], 
                               "lo.y" = smth)
    } else {
      plot_df_lo <- data.frame("lo.x" = x$plot_utilities$reverse_trans(x$plot_utilities$smth$x), 
                               "lo.y" = x$plot_utilities$smth$y)
    }
    
    ggp <- (ggplot(data = plot_df) +
              geom_scattermore(aes_string(x = "m_o", y = "v_o"), alpha = 0.5, 
                               color = "grey30") +
              theme_bw() +
              xlab(paste("Observed (transformed) ", 
                         ifelse(x$plot_utilities$gene_based, "gene average", "gene counts"))) +
              ylab(ifelse(x$plot_utilities$gene_based, "Variance", "log squared-error")) +
              ggtitle("Mean-variance local regression non-parametric fit") +
              geom_line(data = plot_df_lo, aes_string(x = "lo.x", y = "lo.y"), 
                        color = "blue", lwd = 1.4, lty = "solid", alpha = 0.5)
    )
  }else if(x$plot_utilities$method == "none"){
    warning("No plot to be drawn with weighting method = 'none'")
  }else{
    stop("Unknown weighting method")
  }
  
  return(ggp)
  
}
