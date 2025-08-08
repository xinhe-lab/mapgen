
#' @title Make gene Manhattan plot
#'
#' @param gene.pip.res Gene mapping result
#' @param chr Name of the chr column in the gene mapping result
#' @param pos Name of the pos column in the gene mapping result
#' @param gene_name Name of the gene name column in the gene mapping result
#' @param gene_pip Name of the gene PIP column in the gene gene mapping result
#' @param gene.pip.thresh Gene PIP cutoff (default: 0.8)
#' @param highlight Highlight genes with gene PIP > gene.pip.thresh
#' @param ylim Truncate gene PIP to ylim value in the plot.
#' @param point.size Size of the points.
#' @param label.size Size of the labels.
#' @param font.size Font size of the text.
#' @param title Title of the plot
#' @param max.overlaps Exclude text labels that overlap too many things.
#' @import ggplot2
#' @import tidyverse
#' @export
gene_manhattan_plot <- function(gene.pip.res,
                                chr='chr',
                                pos='pos',
                                gene_name='gene_name',
                                gene_pip='gene_pip',
                                gene.pip.thresh = 0.8,
                                highlight = TRUE,
                                ylim = 1.2,
                                point.size = 2,
                                label.size = point.size*2,
                                font.size = 15,
                                max.overlaps = 20,
                                title = '') {

  gene.pip.res <- gene.pip.res %>% dplyr::rename(chr = all_of(chr),
                                                 pos = all_of(pos),
                                                 gene_name = all_of(gene_name),
                                                 gene_pip = all_of(gene_pip))

  # Highlight genes
  gene.pip.res <- gene.pip.res %>%
    dplyr::mutate( is_highlight=(gene_pip >= gene.pip.thresh))

  df <- gene.pip.res %>%

    # Compute chromosome size
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr_len = max(as.numeric(pos))) %>%

    # Calculate cumulative position of each chromosome
    dplyr::mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%

    # Add this info to the initial dataset
    dplyr::left_join(gene.pip.res, ., by=c('chr'='chr')) %>%

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
    scale_color_manual(values = rep(c('grey', 'skyblue'), 22 )) +

    # custom X axis:
    scale_x_continuous(label = gsub('chr','', axis.df$chr, ignore.case = TRUE), breaks= axis.df$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,ylim)) +     # remove space between plot area and x axis

    # Add a dotted line for significant PIP cutoff
    geom_hline(yintercept=gene.pip.thresh, linetype='dashed', color = 'red') +

    # Add highlighted points
    ggrastr::geom_point_rast(data=subset(df, is_highlight==TRUE), color='orange', size=point.size) +

    # Add label using ggrepel to avoid overlapping
    ggrepel::geom_label_repel( data=subset(df, is_highlight==TRUE),
                               aes(label=gene_name),
                               size=label.size,
                               min.segment.length = 0,
                               label.size = NA,
                               fill = alpha(c('white'),0),
                               max.overlaps = max.overlaps,
                               fontface = 'italic') +

    # Custom the theme:
    theme_bw() +
    theme(
      text = element_text(size = font.size),
      legend.position='none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(title) +
    theme() +
    xlab('Chromosome') +
    ylab('Gene PIP')

  return(p)

}

#' @title Make a structure plot of partitioned PIP by locus
#' @description Making a structure plot of partitioned PIP by locus
#' This function is adapted from the 'fastTopics' package
#' https://stephenslab.github.io/fastTopics/
#' @param mat A matrix of the proportion of PIPs partitioned to
#' each annotation category,
#' rows are loci, columns are annotation categories.
#' @param categories annotation categories
#' @param colors Colors of the structure plot categories
#' @param ticks Labels of x-axis
#' @param font.size Font size of structure plot
#' @param highlight Highlight a locus
#' @import ggplot2
#' @export
pip_structure_plot <- function(mat,
                               categories,
                               colors,
                               ticks = NULL,
                               font.size = 9,
                               highlight = NULL,
                               xlab = 'Locus',
                               ylab = 'Proportion',
                               legend.title = 'Category'){

  mat <- na.omit(as.matrix(mat))

  n <- nrow(mat)
  k <- length(categories)
  dat <- data.frame(sample   = rep(1:n,times = k),
                    locus    = rep(rownames(mat),times = k),
                    category = rep(categories,each = n),
                    prop     = c(mat[,categories]))
  dat$category <- factor(dat$category, levels = categories)

  p <- ggplot(dat,aes_string(x = 'sample', y = 'prop',
                             color = 'category',
                             fill = 'category')) +
    geom_col(position = position_fill(reverse = TRUE)) +
    scale_x_continuous(limits = c(0,max(dat$sample) + 1),
                       breaks = ticks,
                       labels = names(ticks)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = xlab, y = ylab, color = legend.title, fill = legend.title) +
    cowplot::theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  if(!is.null(highlight)){
    p <- p + geom_text(aes(label=highlight), y = 1.005,
                       angle = 0, size=3, color = 'black')
  }
  return(p)
}

#' @title Make gene track plot using Gviz
#'
#' @param finemapstats A GRanges object of processed finemapping summary statistics
#' @param region The genomic region to visualize in the format of "chr:start-end".
#' @param gene.annots A GRanges object of gene annotations,
#' needed when plotting loops and gene track.
#' @param txdb A `txdb` object of gene annotations, needed for plotting gene track.
#' @param R LD (correlation) matrix
#' @param LD_snp_ids Variant IDs for the LD matrix.
#' @param bigSNP A `bigsnpr` object attached via bigsnpr::snp_attach()
#' @param counts A list of counts data to display as quantitativ data tracks
#' @param peaks A list of peaks to display as binary data tracks
#' @param loops A list of chromatin loops, e.g. PC-HiC, ABC, etc.
#' @param genome Genome assembly version (default: "hg19").
#' @param filter_loop_genes A vector of gene names. Only show loops connected to the genes.
#' @param filter_loop_snps A vector of SNP IDs. Only show loops connected to the SNPs.
#' @param filter_protein_coding_genes Logical. If TRUE, only shows protein coding gene.
#' @param r2.breaks breaks for LD (r2) levels
#' @param r2.colors colors for LD (r2) levels
#' @param color_piptrack_by Color SNPs in the PIP track by
#' `locus`, `cs` (credible sets), or `none` (same color).
#' @param counts.ylim ylim range (default: between 0 and 1) for the `counts` tracks
#' @param counts.color Colors for the `counts` tracks
#' @param peaks.color Colors for the `peak` tracks
#' @param highlight_snps SNPs (rsIDs) to highlight. Or highlight top SNP ('topSNP').
#' @param highlight.color Colors for the highlighted SNPs
#' @param highlight.width Width for the highlighted locations
#' @param genelabel.side Side to put gene labels,
#' options are: 'above' (default), 'below', 'left', 'right'.
#' @param track.sizes Sizes of the tracks
#' @param rotation.title The rotation angle for the text in the track title
#' @param background.title The background color for the title panel.
#' @param frame If TRUE, plot frames in the panels.
#' @param verbose If TRUE, print detail messages for plotting
#' @import GenomicRanges
#' @import tidyverse
#' @export
track_plot <- function(finemapstats,
                       region,
                       gene.annots = NULL,
                       txdb  = NULL,
                       R = NULL,
                       LD_snp_ids = NULL,
                       bigSNP = NULL,
                       counts  = NULL,
                       peaks  = NULL,
                       loops  = NULL,
                       genome = 'hg19',
                       filter_loop_genes = NULL,
                       filter_loop_snps = NULL,
                       filter_protein_coding_genes = TRUE,
                       r2.breaks = c(0, 0.1, 0.25, 0.75, 0.9, 1),
                       r2.colors = c('black','blue','green','orange','red'),
                       color_piptrack_by = c('locus', 'cs', 'none'),
                       counts.ylim = c(0,1),
                       counts.color,
                       peaks.color = 'navy',
                       loops.color = 'gray',
                       highlight_snps = NULL,
                       highlight.color = 'pink',
                       highlight.width = 1000,
                       genelabel.side = c('above', 'below', 'left', 'right'),
                       track.sizes,
                       rotation.title = 0,
                       background.title = "white",
                       frame = FALSE,
                       verbose = FALSE,
                       ...) {

  genelabel.side <- match.arg(genelabel.side)
  color_piptrack_by <- match.arg(color_piptrack_by)

  seqlevelsStyle(finemapstats) <- 'UCSC'

  if( min(finemapstats$pval) >=0 && max(finemapstats$pval) <= 1 ){
    if(verbose){
      cat('Convert GWAS p-value to -log10(pvalue). \n')
    }
    finemapstats$pval <- -log10(finemapstats$pval)
  }

  # Limit to genomic region to visualize
  region <- as(region, 'GRanges')
  seqlevelsStyle(region) <- 'UCSC'

  rows.in <- which((as.character(seqnames(finemapstats)) == as.character(seqnames(region))) &
                     (finemapstats$pos >= start(region)) &
                     (finemapstats$pos <= end(region)))

  curr.finemapstats <- as.data.frame(finemapstats[rows.in,  ])

  # p-value track
  if (!is.null(R)) {
    # Color SNPs by r2
    if (is.null(LD_snp_ids)) {
      stop("LD_snp_ids is required when using R as input.")
    }

    if (length(r2.breaks) != length(r2.colors)+1) {
      stop("length(r2.breaks) != length(r2.colors) + 1!")
    }

    r2.labels <- sapply(1:(length(r2.breaks)-1), function(x){
      sprintf("%s-%s", r2.breaks[x], r2.breaks[x+1])
    })
    names(r2.colors) <- r2.labels

    curr.finemapstats <- add_LD_from_R(curr.finemapstats, R, LD_snp_ids, r2.breaks, r2.labels)
    if(verbose){ cat(nrow(curr.finemapstats), 'snps included.\n')}

    pval.df <- curr.finemapstats %>%
      dplyr::select(chr, pos, pval, r2.brackets) %>%
      dplyr::mutate(start = pos, end = pos) %>%
      dplyr::select(-pos) %>%
      tidyr::pivot_wider(names_from = r2.brackets, names_sort = TRUE, values_from = 'pval')
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- 'UCSC'
    # r2.groups <- names(mcols(pval.gr))
    r2.groups <- factor(names(mcols(pval.gr)), levels = r2.labels)

    pval.track <- Gviz::DataTrack(range = pval.gr,
                                  genome = genome,
                                  groups = r2.groups,
                                  col = r2.colors,
                                  name = '-log10 p-value')
  } else if (!is.null(bigSNP)) {
    # Color SNPs by r2
    if (length(r2.breaks) != length(r2.colors)+1) {
      stop("length(r2.breaks) != length(r2.colors) + 1!")
    }

    r2.labels <- sapply(1:(length(r2.breaks)-1), function(x){
      sprintf("%s-%s", r2.breaks[x], r2.breaks[x+1])
    })
    names(r2.colors) <- r2.labels

    curr.finemapstats <- add_LD_from_bigSNP(curr.finemapstats, bigSNP, r2.breaks, r2.labels)
    if(verbose){ cat(nrow(curr.finemapstats), 'snps included.\n')}

    pval.df <- curr.finemapstats %>% dplyr::select(chr, pos, pval, r2.brackets) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.df <- pval.df %>% tidyr::pivot_wider(names_from = r2.brackets, names_sort = TRUE, values_from = 'pval')
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- 'UCSC'
    # r2.groups <- names(mcols(pval.gr))
    r2.groups <- factor(names(mcols(pval.gr)), levels = r2.labels)

    pval.track <- Gviz::DataTrack(range = pval.gr,
                                  genome = genome,
                                  groups = r2.groups,
                                  col = r2.colors,
                                  name = '-log10 p-value')
  } else {
    if(verbose){ cat(nrow(curr.finemapstats), 'snps included.\n')}
    pval.df <- curr.finemapstats %>% dplyr::select(chr, pos, pval) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- 'UCSC'

    pval.track <- Gviz::DataTrack(range = pval.gr,
                                  genome = genome,
                                  name = '-log10 p-value',
                                  col = 'black')
  }

  dpars.pval <- list(col.title = 'black',
                     col.axis = 'black',
                     # col.border.title = 'lightgray',
                     # col.frame = 'lightgray',
                     frame = frame,
                     # rotation.title = rotation.title,
                     cex.axis = 0.6)
  Gviz::displayPars(pval.track) <- dpars.pval

  # PIP track
  if(color_piptrack_by == 'locus'){
    if(verbose){ cat('Color SNPs in PIP track by loci. \n')}
    pip.df <- curr.finemapstats %>% dplyr::select(chr, pos, pip, locus) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pip.df <- pip.df %>% tidyr::pivot_wider(names_from = locus, names_sort = TRUE, values_from = 'pip')
    pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
    seqlevelsStyle(pip.gr) <- 'UCSC'
    pip.groups <- names(mcols(pip.gr))
  }else if(color_piptrack_by == 'cs'){
    if(verbose){ cat('Color SNPs in PIP track by credible sets.\n')}
    pip.df <- curr.finemapstats %>% dplyr::select(chr, pos, pip, cs) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pip.df <- pip.df %>% tidyr::pivot_wider(names_from = cs, names_sort = TRUE, values_from = 'pip')
    pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
    seqlevelsStyle(pip.gr) <- 'UCSC'
    pip.groups <- names(mcols(pip.gr))
  }else{
    pip.df <- curr.finemapstats %>% dplyr::select(chr, pos, pip) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
    seqlevelsStyle(pip.gr) <- 'UCSC'
    pip.groups <- NULL
  }

  pip.track <- Gviz::DataTrack(range = pip.gr,
                               genome = genome,
                               groups = pip.groups,
                               name = 'PIP',
                               legend = FALSE)

  dpars.pip <- list(col.title = 'black',
                    col.axis = 'black',
                    # col.border.title = 'lightgray',
                    # col.frame = 'lightgray',
                    frame = frame,
                    # rotation.title = rotation.title,
                    cex.axis = 0.6)
  Gviz::displayPars(pip.track) <- dpars.pip

  # Data tracks
  if (length(counts) > 0) {
    dpars.data <- list(col.title = 'black',
                       col.axis = 'black',
                       # col.border.title = 'lightgray',
                       # col.frame = 'lightgray',
                       frame = frame,
                       rotation.title = rotation.title,
                       cex.axis = 0.2)

    if(missing(counts.color)){
      counts.color <- seq_along(counts)
    }

    counts.tracks <- lapply(names(counts), function(x){
      if(verbose){ cat('Adding', x, 'track...\n') }

      counts.gr <- counts[[x]]
      i <- which(names(counts) == x)
      seqlevelsStyle(counts.gr) <- 'UCSC'
      seqlevels(counts.gr, pruning.mode = 'coarse') <- paste0('chr',1:22)
      counts.track <- Gviz::DataTrack(range = counts.gr,
                                      type = 'h',
                                      genome = genome,
                                      col = counts.color[i],
                                      name = x,
                                      showAxis=FALSE,
                                      ylim = counts.ylim)

      Gviz::displayPars(counts.track) <- dpars.data
      counts.track
    })
  }else{
    counts.tracks <- NULL
  }

  # Peak annotation tracks
  if (length(peaks) > 0) {
    dpars.peaks <- list(col.title = 'black',
                        col.axis = 'black',
                        # col.border.title = 'lightgray',
                        # col.frame = 'lightgray',
                        frame = frame,
                        rotation.title = rotation.title,
                        cex.axis = 0.2)

    if(missing(peaks.color)){
      peaks.color <- seq_along(peaks)
    }

    peaks.tracks <- lapply(names(peaks), function(x){
      if(verbose){ cat('Adding', x, 'track...\n') }
      peaks.gr <- peaks[[x]]
      i <- which(names(peaks) == x)
      seqlevelsStyle(peaks.gr) <- 'UCSC'
      seqlevels(peaks.gr, pruning.mode = 'coarse') <- paste0('chr',1:22)
      peaks.gr$score <- 1
      peaks.gr <- peaks.gr[countOverlaps(peaks.gr, region) > 0]
      peaks.track <- Gviz::DataTrack(range = peaks.gr,
                                     type = 'h',
                                     genome = genome,
                                     name = x,
                                     col = peaks.color[i],
                                     showAxis = FALSE,
                                     ylim = c(0,1))
      Gviz::displayPars(peaks.track) <- dpars.peaks
      peaks.track
    })
  }else{
    peaks.tracks <- NULL
  }

  # Chromatin loop tracks
  if (length(loops) > 0) {
    dpars.loops <- list(col.interactions = loops.color,
                        col.anchors.fill = 'blue',
                        col.anchors.line = 'black',
                        interaction.dimension = 'width',
                        interaction.measure = 'score',
                        plot.trans = FALSE,
                        plot.outside = FALSE,
                        col.outside='lightblue',
                        anchor.height = 0.1,
                        col.title = 'black',
                        col.axis = 'black',
                        # col.border.title = 'lightgray',
                        # col.frame = 'lightgray',
                        frame = frame,
                        rotation.title = rotation.title,
                        cex.axis = 0.2)

    if (is.null(gene.annots)){
      stop("'gene.annots' is needed for plotting loops!")
    }
    gene.annots$chr <- as.character(seqnames(gene.annots))
    gene.annots$tss <- start(resize(gene.annots, width = 1))

    loops.tracks <- lapply(names(loops), function(x){
      if(verbose){ cat('Adding', x, 'track...\n') }
      loops.gr <- loops[[x]]
      if (filter_protein_coding_genes){
        loops.gr <- loops.gr[loops.gr$gene_name %in% gene.annots$gene_name]
      }

      loops_promoters.gr <- GRanges(seqnames = seqnames(loops.gr),
                                    ranges = IRanges(start = loops.gr$promoter_start,
                                                     end = loops.gr$promoter_end),
                                    score = loops.gr$score,
                                    gene = loops.gr$gene_name)
      loops_enhancers.gr <- GRanges(seqnames = seqnames(loops.gr),
                                    ranges = IRanges(start = start(loops.gr),
                                                     end = end(loops.gr)))
      loops.obj <- GenomicInteractions::GenomicInteractions(anchor1 = loops_promoters.gr,
                                                            anchor2 = loops_enhancers.gr)
      loops.obj$counts <- round(loops.obj$anchor1.score)

      if (length(filter_loop_genes) >0 ) {
        cat(sprintf('Only show %s loops linked to gene: %s \n',
                    x, paste(filter_loop_genes, collapse = ',')))
        loops.obj <- loops.obj[which(loops.obj$anchor1.gene %in% filter_loop_genes),]
      }

      if (length(filter_loop_snps) > 0){
        cat(sprintf('Only show %s loops linked to SNP: %s \n',
                    x, paste(filter_loop_snps, collapse = ',')))
        highlighted.snps.gr <- finemapstats[finemapstats$snp %in% filter_loop_snps]
        loops.obj <- subsetByOverlaps(loops.obj, highlighted.snps.gr)
      }

      loops.track <- GenomicInteractions::InteractionTrack(loops.obj, name = x)
      Gviz::displayPars(loops.track) <- dpars.loops
      loops.track
    })
  }else{
    loops.tracks <- NULL
  }

  # gene track
  if (!is.null(txdb)) {
    gene.track <- make_genetrack_obj(region,
                                     txdb,
                                     gene.annots,
                                     genome,
                                     name = 'Gene',
                                     filter_protein_coding_genes = filter_protein_coding_genes)

    dpars.genes <- list(col.title = 'black',
                        col.axis = 'black',
                        # col.border.title = 'lightgray',
                        # col.frame = 'lightgray',
                        frame = FALSE,
                        rotation.title = rotation.title,
                        cex.axis = 0.2)
    displayPars(gene.track) <- dpars.genes

    # Restrict to protein coding genes
    if (isTRUE(filter_protein_coding_genes)){
      gene.track <- gene.track[symbol(gene.track) %in% gene.annots$gene_name]
    }

  } else {
    gene.track <- NULL
  }

  # Genome axis track
  axisTrack <- Gviz::GenomeAxisTrack()

  # List all tracks
  list.of.tracks <- c(pval.track,
                      pip.track,
                      counts.tracks,
                      peaks.tracks,
                      loops.tracks,
                      gene.track,
                      axisTrack)

  if (missing(track.sizes)){
    n.counts.tracks <- length(counts.tracks)
    n.peaks.tracks <- length(peaks.tracks)
    n.loops.tracks <- length(loops.tracks)
    if (is.null(gene.track)){
      n.gene.track <- 0
    }else{
      n.gene.track <- 1
    }
    track.sizes <- c(1,
                     0.6,
                     rep(0.3, n.counts.tracks),
                     rep(0.2, n.peaks.tracks),
                     rep(0.6, n.loops.tracks),
                     rep(0.5, n.gene.track),
                     0.3)
  }

  # Highlight SNPs
  if (length(highlight_snps) > 0){

    if('topSNP' %in% highlight_snps){
      highlight_snps[which(highlight_snps == "topSNP")] <- curr.finemapstats$snp[which.max(curr.finemapstats$pip)]
    }

    highlight_snps <- intersect(highlight_snps, curr.finemapstats$snp)
    highlight.pos <- curr.finemapstats$pos[match(highlight_snps, curr.finemapstats$snp)]
    cat('Highlight SNP:', highlight_snps[which(highlight_snps %in% curr.finemapstats$snp)], '\n')
    cat('Highlight position:', highlight.pos, '\n')
    if (length(highlight.color) > 1){
      highlight.color <- highlight.color[which(highlight_snps %in% curr.finemapstats$snp)]
    }
    list.of.tracks <- Gviz::HighlightTrack(trackList = list.of.tracks,
                                           start = c(highlight.pos-highlight.width/2),
                                           width = highlight.width,
                                           chromosome = as.character(seqnames(region)),
                                           col = 'white',
                                           fill = highlight.color)
  }

  if(verbose){ cat('Making track plot...\n') }

  Gviz::plotTracks(list.of.tracks,
                   chromosome = as.character(seqnames(region)),
                   transcriptAnnotation = 'symbol',
                   collapseTranscripts = 'longest',
                   from = start(region),
                   to = end(region),
                   sizes = track.sizes,
                   just.group = genelabel.side,
                   panel.only = FALSE,
                   background.title = background.title,
                   col.title = "black",
                   col.axis = "black",
                   ...)

}


# Make gene track object for Gviz
make_genetrack_obj <- function(region,
                               txdb = NULL,
                               gene.annots = NULL,
                               genome = 'hg19',
                               keytype = 'ENSEMBL',
                               name = 'Gene',
                               filter_protein_coding_genes = TRUE){

  if(!is.null(txdb)){
    # use supplied txdb if available.
    cat('Making gene track object using txdb database ...\n')
    grtrack <- Gviz::GeneRegionTrack(range = txdb,
                                     genome = genome,
                                     chromosome = as.character(seqnames(region)),
                                     start = start(region),
                                     end = end(region),
                                     name = name)

    symbol(grtrack) <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                             keys=sub('\\.\\d.*$', '', gene(grtrack)),
                                             keytype=keytype,
                                             column='SYMBOL')

    symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), '', symbol(grtrack))

  }else if(!is.null(gene.annots)){
    # use supplied gene.annots if available.
    cat('txdb not available. Making gene track object using gene.annots ...\n')
    grtrack <- Gviz::GeneRegionTrack(range = gene.annots,
                                     genome = genome,
                                     chromosome = as.character(seqnames(region)),
                                     start = start(region),
                                     end = end(region),
                                     name = name,
                                     symbol = gene.annots$gene_name)
  }

  # restrict to protein coding genes
  if (isTRUE(filter_protein_coding_genes)){
    grtrack <- grtrack[symbol(grtrack) %in% gene.annots$gene_name]
  }

  return(grtrack)
}

