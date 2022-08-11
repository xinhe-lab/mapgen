
#' @title Make gene manhattan plot
#'
#' @param gene.pip.res Gene level finemapping result
#' @param chr Name of the chr column in the gene finemapping summary statistics data
#' @param pos Name of the pos column in the gene finemapping summary statistics data
#' @param gene_name Name of the gene name column in the gene finemapping summary statistics data
#' @param gene_pip Name of the gene PIP column in the gene finemapping summary statistics data
#' @param sig.pip Signficant gene PIP (default: 0.8)
#' @param highlight Highlight genes with gene PIP > sig.pip
#' @param ylim Truncate gene PIP to ylim value in the plot (default: 1.25).
#' @param point.size Size of the points.
#' @param label.size Size of the labels.
#' @param font.size Font size of the text.
#' @param title Title of the plot
#' @param max.overlaps Exclude text labels that overlap too many things (default: 10).
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
gene_manhattan_plot <- function(gene.pip.res,
                                chr='chr',
                                pos='pos',
                                gene_name='gene_name',
                                gene_pip='gene_pip',
                                sig.pip = 0.8,
                                highlight = TRUE,
                                ylim = 1.25,
                                point.size = 2,
                                label.size = point.size*2,
                                font.size = 15,
                                max.overlaps = 10,
                                title = '') {

  gene.pip.res <- gene.pip.res %>% dplyr::rename(chr = all_of(chr),
                                                 pos = all_of(pos),
                                                 gene_name = all_of(gene_name),
                                                 gene_pip = all_of(gene_pip))

  # Highlight genes
  gene.pip.res <- gene.pip.res %>%
    dplyr::mutate( is_highlight=(gene_pip >= sig.pip))

  df <- gene.pip.res %>%

    # Compute chromosome size
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr_len = max(as.numeric(pos))) %>%

    # Calculate cumulative position of each chromosome
    dplyr::mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%

    # Add this info to the initial dataset
    dplyr::left_join(gene.pip.res, ., by=c("chr"="chr")) %>%

    # Add a cumulative position of each SNP
    dplyr::arrange(chr, pos) %>%
    dplyr::mutate( pos_cum=pos+tot)

  axis.df <- df %>%
    dplyr::group_by(chr) %>%
    dplyr::summarize(center=( max(pos_cum) + min(pos_cum) ) / 2 )

  p <- ggplot(df,
              aes(x = pos_cum, y = gene_pip)) +

    # Show all points
    ggrastr::geom_point_rast( aes(color=as.factor(chr)), size=point.size) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +

    # custom X axis:
    scale_x_continuous(label = gsub("chr","", axis.df$chr, ignore.case = TRUE), breaks= axis.df$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,ylim)) +     # remove space between plot area and x axis

    # Add highlighted points
    ggrastr::geom_point_rast(data=subset(df, is_highlight==TRUE), color="orange", size=point.size) +

    # Add label using ggrepel to avoid overlapping
    ggrepel::geom_label_repel( data=subset(df, is_highlight==TRUE),
                               aes(label=gene_name),
                               size=label.size,
                               min.segment.length = 0,
                               label.size = NA,
                               fill = alpha(c("white"),0),
                               max.overlaps = max.overlaps) +

    # Custom the theme:
    theme_bw() +
    theme(
      text = element_text(size = font.size),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(title) +
    theme() +
    xlab("Chromosome") +
    ylab("Gene PIP")

  return(p)

}



#' @title structure plot
#' @description this function is adapted from the 'fastTopics' package
#' https://stephenslab.github.io/fastTopics/
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
  p <- ggplot(dat,aes_string(x = "sample",y = "prop",
                             color = "category",
                             fill = "category")) +
    geom_col() +
    scale_x_continuous(limits = c(0,max(dat$sample) + 1),
                       breaks = ticks,
                       labels = names(ticks)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "",y = "Proportion") +
    cowplot::theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1))

  if(!is.null(highlight)){
    p <- p + geom_text(aes(label=highlight), y = 1.005, angle = 0, size=3, color = "black")
  }
  return(p)
}

#' @title Compile input data for making structure plot
#' @description this function is adapted from the 'fastTopics' package
#' https://stephenslab.github.io/fastTopics/
#'
#' @param mat matrix of input data
#' @param categories annotation categories
#' @return a data frame for making structure plot
#' @export
#'
compile_structure_plot_data <- function (mat, categories) {
  mat <- na.omit(as.matrix(mat))
  n <- nrow(mat)
  k <- length(categories)
  dat <- data.frame(sample   = rep(1:n,times = k),
                    locus    = rep(rownames(mat),times = k),
                    category = rep(categories,each = n),
                    prop     = c(mat[,categories]))
  dat$category <- factor(dat$category, levels = categories)
  return(dat)
}
