#'Spaghetti plot for Specific Gene Set
#'
#' @param gs_index index of the specific gene set in \code{gmt}. 
#' @param gmt a \code{list} of elements : geneset, geneset.name and geneset.description.
#' @param expr_mat a data frame with numerics of size \code{G x n} contraining the raw RNA-seq counts from \code{n} samples for \code{G} genes.
#' @param design a data frame containing the information of each sample (SampleID).
#' @param var_time the \code{time} or \code{visit} variable contained in \code{design} data.
#' @param var_patient the patient variable contained in \code{design} data.
#' @param var_group a group variable in \code{design} data to divide into two facets. Default is \code{NULL}.
#' @param var_subgroup a subgroup variable in \code{design} data to add 2 curves on plot for each subgroup. Default is \code{NULL}.
#' @param plotChoice to choose which type of plot. Default is \code{c("Patient medians", "Individual")}. 
#' 
#' @return 
#' @import dplyr
#' @import ggplot2
#' @import rlang
#' @export
#' 
#' @examples
#' data(baduel_5gs) 
#' 
#' expr_norm_corr <- as.data.frame(expr_norm_corr)
#' #Remove "Vernalized" columns 
#' expr_norm_corr <- expr_norm_corr %>% select(!ends_with("V")) 
#' 
#' #Change name of sample 
#' colnames(design)[1] <- "SampleID"
#' design <- design %>% filter(Vernalized == FALSE)
#' #Modify the Population ID to obtain value of population for each replicate
#' design <- design %>% mutate(PopulationID = Population:Replicate)
#' 
#' dearseq::spaghettiPlot_1GS(gs_index = 3, gmt = baduel_gmt, expr_mat = expr_norm_corr, design = design, var_time = AgeWeeks, var_patient = PopulationID)

spaghettiPlot_1GS <- function(gs_index, gmt, expr_mat, design, var_time, var_patient, var_group = NULL, var_subgroup = NULL, plotChoice = c("Patient medians", "Individual")){ 

  #Before to call a tidy evaluation function (ggplot) inside of another function
  # use enquo() and !! before object of the function
  var_subgroup_tidy <- enquo(var_subgroup)
  var_patient_tidy <- enquo(var_patient)
  var_time_tidy <- enquo(var_time)
  var_group_tidy <- enquo(var_group)

  #Select expression of genes in genesets
  gs_expr <- expr_mat %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    filter(SYMBOL %in% gmt$genesets[[gs_index]]) 

  data2plot <- reshape2::melt(gs_expr, value.name = "Normalized_Expression", variable.name = "SampleID") %>% 
    left_join(y = design) %>% 
    mutate(SYMBOL = as.factor(SYMBOL))

  data2plot <- data2plot %>% mutate(group_line = SYMBOL:!!var_patient_tidy)
  
  data2plot <- data2plot %>% 
    group_by(SYMBOL, !!var_group_tidy, !!var_subgroup_tidy, !!var_patient_tidy) %>% 
    mutate(Standardized_Expression = (Normalized_Expression - mean(Normalized_Expression, na.rm = TRUE))/sd(Normalized_Expression))

  data2plot_summ <- data2plot %>% ungroup() %>% 
    group_by(SYMBOL, !!var_group_tidy, !!var_subgroup_tidy, !!var_time_tidy) %>%
    summarise(Standardized_Expression = median(Standardized_Expression), 
              .groups = "drop_last") #%>% 

  # to have the vector of time
  var_time_vec <- eval_tidy(var_time_tidy, design)
  
  #Plot of Patient medians
  p1 <- ggplot(data2plot_summ, aes(x = !!var_time_tidy, y = Standardized_Expression)) +
    geom_line(aes(color = SYMBOL, linetype = !!var_subgroup_tidy)) +
    geom_smooth(aes(linetype = !!var_subgroup_tidy), color="black", se=FALSE) +
    theme_bw() +
    scale_x_continuous(breaks = unique(var_time_vec)) +
    scale_color_discrete("Gene") + 
    guides(linetype = guide_legend(override.aes = list(size=1))) +
    ylab("Standardized expression") +
    ggtitle("Patient medians")


  #Plot of Individual data
  p2 <- ggplot(data2plot, aes(x = !!var_time_tidy, y = Standardized_Expression)) +
    geom_line(aes(group = group_line, color = SYMBOL, linetype = !!var_subgroup_tidy), alpha = 0.4) +
    geom_smooth(formula="y~x", method="loess", aes(linetype = !!var_subgroup_tidy), color="black", se=FALSE) +
    theme_bw() +
    scale_x_continuous(breaks = unique(var_time_vec)) +
    ylab("Standardized expression") +
    ggtitle("Individual data")
  
  #If there is var_group, add facet_wrap
  if(!is.null(eval_tidy(var_group_tidy, data = design))){
    p1 <- p1 + facet_wrap(var_group_tidy) 
    p2 <- p2 + facet_wrap(var_group_tidy) 
  }

  #Combine the 2 plots 
  pall <- p1/(p2 + guides(color = "none", linetype = "none")) + plot_layout(guides = 'collect') + 
    plot_annotation(title = paste0(gmt$geneset.name[[gs_index]], ": ", gmt$geneset.description[[gs_index]]))
  
  #Return plot according to plotChoice
  if(identical(plotChoice, "Patient medians")){
    p1 <- p1 +
          plot_annotation(title = paste0(gmt$geneset.name[[gs_index]], ": ", gmt$geneset.description[[gs_index]]))
    return(p1)
  }else  if(identical(plotChoice, "Individual")){
    p2 <- p2 +    
          scale_color_discrete("Gene") + 
          guides(linetype = guide_legend(override.aes = list(size=1))) +
          plot_annotation(title = paste0(gmt$geneset.name[[gs_index]], ": ", gmt$geneset.description[[gs_index]]))
    return(p2)
  }else if(length(plotChoice) == 2){
    return(pall)
  }else{
    warning("You must choose the type of plot to return.")
  }

}

