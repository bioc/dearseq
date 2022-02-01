#' Spaghetti plot for Specific Gene Set
#'
#' @param gs_index index of the specific gene set in \code{gmt}. 
#' @param gmt a \code{list} of elements: \code{geneset}, \code{geneset.name} and
#'  \code{geneset.description} (see \code{\link[GSA]{GSA.read.gmt}}).
#' @param expr_mat a \code{data.frame} with numerics of size \code{G x n} 
#' contraining the raw RNA-seq counts from \code{n} samples for \code{G} genes.
#' @param design a \code{data.frame} or \code{DFrame} containing the information 
#' of each sample (SampleID).
#' @param var_time the \code{time} or \code{visit} variable contained in 
#' \code{design}.
#' @param var_indiv the patient variable contained in \code{design} data.
#' @param sampleIdColname a character string indicating the name of the sample ID 
#' variable in \code{design} to be matched with the \code{colnames} of \code{expr_mat}
#' @param var_group a group variable in \code{design} data to divide into two 
#' facets. Default is \code{NULL}.
#' @param var_subgroup a subgroup variable in \code{design} data to add 2 curves 
#' on plot for each subgroup. Default is \code{NULL}.
#' @param plotChoice to choose which type of plot (either \code{"Medians"}, 
#' \code{"Individual"} or both). Default is \code{c("Medians", "Individual")}. 
#' @param loess_span smoothing span. Default is \code{0.75}.
#' 
#' @return a ggplot2 plot object 
#'
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @importFrom dplyr filter left_join
#' @importFrom rlang eval_tidy
#' @importFrom stats median sd
#' @import ggplot2
#' @export
#' 
#' @examples
#' data(baduel_5gs) 
#' design$Indiv <- design$Population:design$Replicate
#' design$Vern <- ifelse(design$Vernalized, "Vernalized", "Non-vernalized")
#' 
#' spaghettiPlot1GS(gs_index = 3, gmt = baduel_gmt, expr_mat = log2(expr_norm_corr+1), 
#'   design = design, var_time = AgeWeeks, var_indiv = Indiv, 
#'   sampleIdColname = "sample", var_group=Vern, var_subgroup=Population, 
#'   plotChoice = "Medians", loess_span= 1.5) +
#'   xlab("Age (weeks)") + guides(color= "none")

spaghettiPlot1GS <- function(gs_index, gmt, expr_mat, design, var_time, 
                             var_indiv, sampleIdColname, 
                             var_group = NULL, var_subgroup = NULL, 
                             plotChoice = c("Medians", "Individual"),
                             loess_span = 0.75){ 

  if(is.matrix(expr_mat)){
    expr_mat <- as.data.frame(expr_mat)
  }
  stopifnot(is.data.frame(expr_mat))
  if(class(design) == "DFrame"){
    design <- as.data.frame(design)
  }
  stopifnot(is.data.frame(design))
  
  #Before to call a tidy evaluation function (ggplot) inside of another function
  # use enquo() and !! before object of the function
  var_subgroup_tidy <- enquo(var_subgroup)
  var_indiv_tidy <- enquo(var_indiv)
  var_time_tidy <- enquo(var_time)
  var_group_tidy <- enquo(var_group)
  
  SYMBOL_tidy <- sym("SYMBOL")
  Normalized_Expression_tidy <- sym("Normalized_Expression")
  Standardized_Expression_tidy <- sym("Standardized_Expression")

  #Select expression of genes in genesets
  gs_expr <- expr_mat %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    dplyr::filter(!!SYMBOL_tidy %in% gmt$genesets[[gs_index]]) 

  colnames(design)[which(colnames(design) == sampleIdColname)] <- "SampleID"
  data2plot <- reshape2::melt(gs_expr, value.name = "Normalized_Expression", 
                              variable.name = "SampleID", id.vars = "SYMBOL") %>% 
    dplyr::left_join(y = design, by = "SampleID") %>% 
    dplyr::mutate("SYMBOL" = as.factor(!!SYMBOL_tidy))

  data2plot <- data2plot %>% 
    dplyr::mutate("group_line" = !!SYMBOL_tidy:!!var_indiv_tidy)
  
  data2plot <- data2plot %>% 
    dplyr::group_by(!!SYMBOL_tidy, !!var_group_tidy, !!var_subgroup_tidy, 
                    !!var_indiv_tidy) %>% 
    dplyr::mutate("Standardized_Expression" = ((!!Normalized_Expression_tidy - 
                                               mean(!!Normalized_Expression_tidy, 
                                                    na.rm = TRUE)) / 
                                               sd(!!Normalized_Expression_tidy)))

  data2plot_summ <- data2plot %>% dplyr::ungroup() %>% 
    dplyr::group_by(!!SYMBOL_tidy, !!var_group_tidy, !!var_subgroup_tidy, !!var_time_tidy) %>%
    dplyr::summarise("Standardized_Expression" = median(!!Standardized_Expression_tidy), 
              .groups = "drop_last") #%>% 

  # to have the vector of time
  var_time_vec <- rlang::eval_tidy(var_time_tidy, design)
  
  #Plot of Patient medians
  p1 <- ggplot(data2plot_summ, aes_string(x = var_time_tidy, y = "Standardized_Expression")) +
    geom_line(aes_string(color = "SYMBOL", linetype = var_subgroup_tidy)) +
    geom_smooth(aes_string(linetype = var_subgroup_tidy), color="black", 
                se=FALSE, method = 'loess', formula = 'y ~ x', span = loess_span) +
    theme_bw() +
    scale_x_continuous(breaks = unique(var_time_vec)) +
    scale_color_discrete("Gene") + 
    guides(linetype = guide_legend(override.aes = list(size=1))) +
    ylab("Standardized expression") +
    ggtitle("Median across individuals")


  #Plot of Individual data
  p2 <- ggplot(data2plot, aes_string(x = var_time_tidy, 
                                     y = "Standardized_Expression")) +
    geom_line(aes_string(group = "group_line", color = "SYMBOL", 
                         linetype = var_subgroup_tidy), alpha = 0.4) +
    geom_smooth(aes_string(linetype = var_subgroup_tidy), color="black", 
                se=FALSE, formula="y~x", method="loess", span = loess_span) +
    theme_bw() +
    scale_x_continuous(breaks = unique(var_time_vec)) +
    ylab("Standardized gene expression") +
    ggtitle("Individual data")
  
  #If there is var_group, add facet_wrap
  if(!is.null(eval_tidy(var_group_tidy, data = design))){
    p1 <- p1 + facet_wrap(var_group_tidy) 
    p2 <- p2 + facet_wrap(var_group_tidy) 
  }

  #Combine the 2 plots 
  pall <- p1/(p2 + guides(color = "none", linetype = "none")) + 
    plot_layout(guides = 'collect') + 
    plot_annotation(title = paste0(gmt$geneset.name[[gs_index]], ": ", 
                                   gmt$geneset.description[[gs_index]]))
  
  #Return plot according to plotChoice
  if(identical(plotChoice, "Medians")){
    p1 <- p1 +
          plot_annotation(title = paste0(gmt$geneset.name[[gs_index]], ": ", 
                                         gmt$geneset.description[[gs_index]]))
    return(p1)
  }else  if(identical(plotChoice, "Individual")){
    p2 <- p2 +    
          scale_color_discrete("Gene") + 
          guides(linetype = guide_legend(override.aes = list(size=1))) +
          plot_annotation(title = paste0(gmt$geneset.name[[gs_index]], ": ", 
                                         gmt$geneset.description[[gs_index]]))
    return(p2)
  }else if(length(plotChoice) == 2){
    return(pall)
  }else{
    warning("You must choose the type of plot to return.")
  }

}

