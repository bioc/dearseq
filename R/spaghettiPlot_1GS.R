#'Spaghetti plot for Specific Gene Set
#'
#' @param gs_index index of the gene set in the \code{gmt} 
#' @param expr_mat where rownames = names of genes and colames = names of sampleId
#' @param gmt list of element whose gene sets, their name and description
#' @param design sampleinfo
#' @param var_time in \code{design}
#' @param var_patient in \code{design}
#' @param var_group in \code{design} (Vaccine)
#' @param var_add in \code{design} and in example of Covax is Recovered
#' 
#' @return 
#' @export
#' 
#' @examples

spaghettiPlot_1GS <- function(gs_index, expr_mat, gmt, design, var_time, 
                              var_patient, var_group, var_add){ 
  
  gs_expr <- expr_mat %>% 
    tibble::rownames_to_column(var="SYMBOL") %>% 
    filter(SYMBOL %in% gmt$genesets[[gs_index]]) 
  
  data2plot <- reshape2::melt(gs_expr, value.name = "Normalized_Expression", variable.name = "SampleID") %>% 
    left_join(y = design) %>% 
    mutate(SYMBOL = as.factor(SYMBOL))
  data2plot$group_line <- data2plot$SYMBOL:data2plot$Patient
  
  
  data2plot$Recovered <- as.factor(data2plot$Recovered)
  levels(data2plot$Recovered) <- c(paste0("No ", "n=", recovered_no_n),
                                   paste0("Yes ", "n=", recovered_yes_n))
  levels(data2plot$Vaccine) <- c(paste0("Pfizer ", "n=", pfizer_n),
                                 paste0("Moderna ", "n=", moderna_n))
  
  
  data2plot <- data2plot %>% 
    group_by(SYMBOL, Vaccine, Recovered, Patient) %>% 
    mutate(Standardized_Expression = (Normalized_Expression - mean(Normalized_Expression, na.rm = TRUE))/sd(Normalized_Expression))
  
  data2plot_summ <- data2plot %>% ungroup() %>% 
    group_by(SYMBOL, Vaccine, Recovered, Time) %>%
    summarise(Standardized_Expression=median(Standardized_Expression), 
              .groups = "drop_last") #%>% 
  #mutate(Standardized_expression = (Median_Normalized_Expression - mean(Median_Normalized_Expression))/sd(Median_Normalized_Expression))
  p1 <- ggplot(data2plot_summ, aes(x=Time, y = Standardized_Expression)) +
    geom_line(aes(color = SYMBOL, linetype = Recovered)) + 
    geom_smooth(aes(linetype = Recovered), color="black", se=FALSE) +
    facet_wrap(~Vaccine) +
    theme_bw() +
    scale_x_continuous(breaks=unique(design$Time)) +
    scale_color_discrete("Gene") + 
    guides(linetype=guide_legend(override.aes = list(size=1))) +
    ylab("Standardized expression") +
    ggtitle("Patient medians")
  
  p2 <- ggplot(data2plot, aes(x=Time, y = Standardized_Expression)) +
    geom_line(aes(group = group_line, color = SYMBOL, linetype = Recovered), alpha=0.4) + 
    geom_smooth(formula="y~x", method="loess", aes(linetype = Recovered), color="black", se=FALSE) +
    facet_wrap(~Vaccine) +
    theme_bw() +
    scale_x_continuous(breaks=unique(design$Time)) +
    ylab("Standardized expression") +
    ggtitle("Individual data")
  
  pall <- p1/(p2 + guides(color="none", linetype="none")) + plot_layout(guides = 'collect') + 
    plot_annotation(title=paste0(gmt$geneset.names[[gs_index]], ": ", gmt$genesets.descriptions[[gs_index]]))
  return(pall)