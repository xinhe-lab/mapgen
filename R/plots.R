#' @title Make gene manhattan plot
#'
#' @param gene.pip Gene level finemapping result
#' @param chr Name of the chr column in the gene finemapping summary statistics data
#' @param POS Name of the pos column in the gene finemapping summary statistics data
#' @param gene Name of the gene name column in the gene finemapping summary statistics data
#' @param pip Name of the gene PIP column in the gene finemapping summary statistics data
#' @param sig.pip Signficance cutoff of gene PIP
#' @param annotate.pip annotate genes with gene PIP > annotate.pip
#' @param highlight.genes Genes to highlight
#' @param ylim Truncate gene PIP to ylim value in the plot
#' @param title Title of the plot
#' @param max.overlaps Limit of overlapping labels (ggrepel setting)
#' @import ggplot2
#' @import tidyverse
#' @export
#'
gene_manhattan_plot <- function(gene.pip,
                                chr='chr',
                                pos='pos',
                                gene='gene_name',
                                pip='gene_pip',
                                highlight = FALSE,
                                annotate = FALSE,
                                sig.pip = 0.5,
                                annotate.pip = NULL,
                                highlight.genes = NULL,
                                ylim = NULL,
                                title = '',
                                max.overlaps = 20) {

  gene.pip <- gene.pip %>% dplyr::rename(chr = all_of(chr),
                                         pos = all_of(pos),
                                         gene = all_of(gene),
                                         pip = all_of(pip))

  # Compute chromosome size and calculate cumulative position of each chromosome
  data_cum <- gene.pip %>%
    group_by(chr) %>%
    summarise(max_pos = max(as.numeric(pos))) %>%
    mutate(pos_add = lag(cumsum(max_pos), default = 0)) %>%
    dplyr::select(chr, pos_add)

  gene.pip <- gene.pip %>%
    inner_join(data_cum, by = 'chr') %>%
    mutate(pos_cum = pos + pos_add)

  # Highlight genes
  if(!is.null(highlight.genes)) {
    gene.pip <- gene.pip %>%
      mutate( is_highlight=ifelse(gene %in% highlight.genes, 'yes', 'no'))
  }else{
    gene.pip$is_highlight <- 'no'
  }

  if(!is.null(annotate.pip)){
    gene.pip <- gene.pip %>%
      mutate( is_annotate=ifelse(pip >= annotate.pip, 'yes', 'no'))
  }else{
    gene.pip$is_annotate <- 'no'
  }

  axis_set <- gene.pip %>%
    group_by(chr) %>%
    summarize(center = mean(pos_cum))

  if(is.null(ylim)){
    ylim <- ceiling(max(gene.pip$pip))
  }

  gene.pip$pip[gene.pip$pip > ylim] <- ylim

  p <- ggplot(gene.pip,
              aes(x = pos_cum, y = pip, color = as_factor(chr), size = pip)) +

    # Show all points
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = sig.pip, color = 'grey40', linetype = 'dashed') +
    scale_color_manual(values = rep(c('#276FBF', '#183059'), unique(length(axis_set$chr)))) +
    scale_size_continuous(range = c(0.5, 1)) +

    # Custom X, Y axis:
    scale_x_continuous(label = gsub('chr','', axis_set$chr, ignore.case = TRUE), breaks = axis_set$center) +
    scale_y_continuous(limits = c(0, ylim)) +

    # Add highlighted points
    geom_point(data=subset(gene.pip, is_highlight=='yes'), color='orange', size=1) +

    # Add label using ggrepel
    geom_label_repel( data=subset(gene.pip, is_annotate=='yes'),
                      aes(label=gene), size=2, max.overlaps = max.overlaps) +

    labs(x = 'Chr',
         y = 'Gene PIP',
         title = title) +

    # Custom the theme:
    theme_minimal() +
    theme(
      legend.position = 'none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 8, vjust = 5)
    )

  return(p)

}


#' @title structure plot
#'
#' @param dat a data frame obtained from @rdname compile_structure_plot_data
#' @param colors Colors of the structure plot categories
#' @param ticks Labels of x-axis
#' @param font.size Font size of structure plot
#' @param highlight Highlight a sample
#' @import ggplot2
#' @export
#'
structure_plot <- function (dat, colors, ticks = NULL,
                            font.size = 9, highlight = NULL){
  p <- ggplot(dat,aes_string(x = "sample",y = "prop",color = "category",
                             fill = "category")) +
    geom_col() +
    scale_x_continuous(limits = c(0,max(dat$sample) + 1),
                       breaks = ticks,
                       labels = names(ticks)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "",y = "Proportion") +
    theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1))

  if(!is.na(highlight)){
    p <- p + geom_text(aes(label=highlight), y = 1.005, angle = 0, size=3, color = "black")
  }
  return(p)
}

# compile input data for making structure plot
compile_structure_plot_data <- function (mat, categories) {
  mat <- na.omit(as.matrix(mat))
  n <- nrow(mat)
  k <- length(categories)
  dat <- data.frame(sample   = rep(1:n,times = k),
                    locus    = rep(rownames(mat),times = k),
                    category = rep(categories,each = n),
                    prop     = c(mat[,categories]))
  dat$category <- factor(dat$category, categories)
  return(dat)
}
