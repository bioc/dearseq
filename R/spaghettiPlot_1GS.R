#'Spaghetti plot for Specific Gene Set
#'
#' @param gs_index index of the gene set in \code{gmt} 
#' @param expr_mat where rownames = names of genes and colnames = names of sampleId
#' @param gmt list of element whose gene sets, their name and description
#' @param design sampleinfo
#' @param var_time in \code{design}
#' @param var_patient in \code{design}
#' @param var_group in \code{design} (Vaccine)
#' @param var_add in \code{design} and in example of Covax is Recovered
#' 
#' @return 
#' @import dplyr
#' @import ggplot2
#' @export
#' 
#' @examples
#' 
#' 

spaghettiPlot_1GS <- function(gs_index, expr_mat, gmt, design, var_time, var_patient, var_group, var_add){ 
  
  #Before to call a tidy evaluation function (ggplot) inside of another function
  # use enquo() and !! before object of the function
  var_add_tidy <- enquo(var_add)
  var_patient_tidy <- enquo(var_patient)
  var_time_tidy <- enquo(var_time)
  var_group_tidy <- enquo(var_group)
  
  gs_expr <- expr_mat %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    filter(SYMBOL %in% gmt$genesets[[gs_index]]) 
  
  data2plot <- reshape2::melt(gs_expr, value.name = "Normalized_Expression", variable.name = "SampleID") %>% 
    left_join(y = design) %>% 
    mutate(SYMBOL = as.factor(SYMBOL))

  data2plot <- data2plot %>% mutate(group_line = SYMBOL:!!var_patient_tidy)
  
  
  # data2plot$Recovered <- as.factor(data2plot$Recovered)
  # levels(data2plot$Recovered) <- c(paste0("No ", "n=", recovered_no_n),
  #                                  paste0("Yes ", "n=", recovered_yes_n))
  # levels(data2plot$Vaccine) <- c(paste0("Pfizer ", "n=", pfizer_n),
  #                                paste0("Moderna ", "n=", moderna_n))
  #
  
  data2plot <- data2plot %>% 
    group_by(SYMBOL, !!var_group_tidy, !!var_add_tidy, !!var_patient_tidy) %>% 
    mutate(Standardized_Expression = (Normalized_Expression - mean(Normalized_Expression, na.rm = TRUE))/sd(Normalized_Expression))
  

  data2plot_summ <- data2plot %>% ungroup() %>% 
    group_by(SYMBOL, !!var_group_tidy, !!var_add_tidy, !!var_time_tidy) %>%
    summarise(Standardized_Expression = median(Standardized_Expression), 
              .groups = "drop_last") #%>% 
  #mutate(Standardized_expression = (Median_Normalized_Expression - mean(Median_Normalized_Expression))/sd(Median_Normalized_Expression))
  
  # to have the vector of time
  var_time_vec <- eval_tidy(var_time_tidy, design)

  p1 <- ggplot(data2plot_summ, aes(x = !!var_time_tidy, y = Standardized_Expression)) +
    geom_line(aes(color = SYMBOL, linetype = !!var_add_tidy)) +
    geom_smooth(aes(linetype = !!var_add_tidy), color="black", se=FALSE) +
    facet_wrap(var_group_tidy) +
    theme_bw() +
    scale_x_continuous(breaks = unique(var_time_vec)) +
    scale_color_discrete("Gene") + 
    guides(linetype = guide_legend(override.aes = list(size=1))) +
    ylab("Standardized expression") +
    ggtitle("Patient medians")
  

  p2 <- ggplot(data2plot, aes(x = !!var_time_tidy, y = Standardized_Expression)) +
    geom_line(aes(group = group_line, color = SYMBOL, linetype = !!var_add_tidy), alpha=0.4) +
    geom_smooth(formula="y~x", method="loess", aes(linetype = !!var_add_tidy), color="black", se=FALSE) +
    facet_wrap(var_group_tidy) +
    theme_bw() +
    scale_x_continuous(breaks = unique(var_time_vec)) +
    ylab("Standardized expression") +
    ggtitle("Individual data")
  
  
  pall <- p1/(p2 + guides(color = "none", linetype = "none")) + plot_layout(guides = 'collect') + 
    plot_annotation(title = paste0(gmt$geneset.names[[gs_index]], ": ", gmt$genesets.descriptions[[gs_index]]))
  return(pall)
}