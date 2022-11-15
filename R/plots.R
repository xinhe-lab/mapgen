
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


#' @title Make gene track plot
#'
#' @param finemapstats A GRanges object or data frame of finemapping summary statistics
#' @param region A GRanges object or data frame for the genomic range to plot
#' @param gene.annots A GRanges object of gene annotations
#' @param bigSNP A bigsnpr object attached via bigsnpr::snp_attach()
#' @param txdb A txdb object of gene annotations
#' @param genome Genome assembly version, hg19 (default) or hg38.
#' @param genetrack_db Select a gene annotation database to use. Options:
#' `txdb`: use the `txdb` objec.
#' `gene.annots`: use the `gene.annots` object.
#' `UCSC` uses `UCSC knownGene` annotations.
#' @param filter_protein_coding_genes If TRUE, only shows protein coding gene
#' @param countsdata A list of counts data
#' @param peaks A list of peaks
#' @param HiC_loops A list of HiC loops, e.g. PC-HiC, ABC, etc.
#' @param filter_HiCloops_genes If TRUE, only shows HiC loops connected to
#' the gene(s)
#' @param filter_HiCloops_snps If TRUE, only shows HiC loops connected to
#' the SNP(s)
#' @param data_colors Colors for the `countsdata` tracks
#' @param data_ylim ylim range for the `countsdata` tracks
#' @param highlight_snps SNPs (rsIDs) to highlight
#' @param highlight_colors Colors for the highlighted SNPs
#' @param genelabel_side Side to put gene labels,
#' options are: above (default), right, left, below
#' @param rotation.title Rotation of the track titles
#' @param verbose if TRUE, print detail messages for plotting
#'
#' @import Gviz
#' @importFrom magrittr %>%
#' @export
finemapping_annot_trackplot <- function(finemapstats,
                                        region,
                                        gene.annots,
                                        bigSNP,
                                        txdb,
                                        genome = c("hg19", "hg38"),
                                        genetrack_db = c("txdb", "gene.annots","UCSC"),
                                        filter_protein_coding_genes = TRUE,
                                        countsdata,
                                        peaks,
                                        HiC_loops,
                                        filter_HiCloops_genes = NULL,
                                        filter_HiCloops_snps = NULL,
                                        data_colors = seq_along(countsdata),
                                        data_ylim = c(0,1),
                                        highlight_snps = NULL,
                                        highlight_colors = "pink",
                                        genelabel_side = c("above", "right", "left", "below"),
                                        rotation.title = 90,
                                        verbose = FALSE) {

  genome <- match.arg(genome)
  genetrack_db <- match.arg(genetrack_db)
  genelabel_side <- match.arg(genelabel_side)

  if(verbose){ cat("Making trackplots ...\n") }

  # Prepare GWAS summary stats
  if( min(finemapstats$pval) >=0 && max(finemapstats$pval) <= 1 ){
    if(verbose){
      cat("Convert GWAS p-value to -log10(pvalue). \n")
    }
    finemapstats$pval <- -log10(finemapstats$pval)
  }
  finemapstats <- as(finemapstats, "GRanges")
  seqlevelsStyle(finemapstats) <- "UCSC"

  # Limit to genomic region to visualize
  region <- as(region, "GRanges")
  seqlevelsStyle(region) <- "UCSC"
  curr_finemapstats <- finemapstats[(as.character(seqnames(finemapstats)) == as.character(seqnames(region))) &
                                      (finemapstats$pos >= start(region)) &
                                      (finemapstats$pos <= end(region)),  ]
  curr_finemapstats <- as.data.frame(curr_finemapstats)

  # p-value track
  # Add LD information
  if(!missing(bigSNP)){
    curr_finemapstats <- add_LD_bigSNP(curr_finemapstats, bigSNP)
    if(verbose){ cat(nrow(curr_finemapstats), "snps included.\n")}

    pval.df <- curr_finemapstats %>% dplyr::select(chr, pos, pval, r2) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.df <- pval.df %>% tidyr::pivot_wider(names_from = r2, names_sort = TRUE, values_from = "pval")
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- "UCSC"
    avail_ld_groups <- names(mcols(pval.gr))

    ld_colors <- c("black","blue","green","orange","red")
    names(ld_colors) <- c("0-0.1","0.1-0.25","0.25-0.75","0.75-0.9","0.9-1")
    pval.track <- DataTrack(range = pval.gr,
                            genome = genome,
                            groups = avail_ld_groups,
                            col = ld_colors[avail_ld_groups],
                            name = "-log10 P",
                            legend = TRUE,
                            box.legend = TRUE)
  }else{
    if(verbose){ cat(nrow(curr_finemapstats), "snps included.\n")}
    pval.df <- curr_finemapstats %>% dplyr::select(chr, pos, pval) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- "UCSC"

    pval.track <- DataTrack(range = pval.gr,
                            genome = genome,
                            name = "-log10 P",
                            col = "black",
                            legend = FALSE)
  }

  dpars.pval <- list(col.title = "black",
                     col.axis = "black",
                     col.border.title = "lightgray",
                     col.frame = "lightgray",
                     # rotation.title = rotation.title,
                     frame = TRUE,
                     cex.axis = 0.6)
  displayPars(pval.track) <- dpars.pval

  # PIP track
  pip.df <- curr_finemapstats %>% dplyr::select(chr, pos, pip, locus) %>%
    dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
  pip.df <- pip.df %>% tidyr::pivot_wider(names_from = locus, names_sort = TRUE, values_from = "pip")
  pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
  seqlevelsStyle(pip.gr) <- "UCSC"
  avail_locus_groups <- names(mcols(pip.gr))

  pip.track <- DataTrack(range = pip.gr,
                         genome = genome,
                         groups = avail_locus_groups,
                         name = "PIP",
                         legend = FALSE)

  dpars.pip <- list(col.title = "black",
                    col.axis = "black",
                    col.border.title = "lightgray",
                    col.frame = "lightgray",
                    # rotation.title = rotation.title,
                    frame = TRUE,
                    cex.axis = 0.6)
  displayPars(pip.track) <- dpars.pip

  # Data tracks
  if(!missing(countsdata) && (length(countsdata) > 0)){
    dpars.data <- list(col.title = "black",
                       col.axis = "black",
                       col.border.title = "lightgray",
                       col.frame = "lightgray",
                       rotation.title = rotation.title,
                       frame = TRUE,
                       cex.axis = 0.2)

    countsdata.tracks <- lapply(names(countsdata), function(x){
      if(verbose){ cat("Adding", x, "track...\n") }

      countsdata.gr <- countsdata[[x]]
      i <- which(names(countsdata) == x)
      seqlevelsStyle(countsdata.gr) <- "UCSC"
      seqlevels(countsdata.gr, pruning.mode = "coarse") <- paste0("chr",1:22)
      countsdata.track <- DataTrack(range = countsdata.gr,
                                    type = 'h',
                                    genome = genome,
                                    col = data_colors[i],
                                    name = x,
                                    showAxis=FALSE,
                                    ylim = data_ylim)

      displayPars(countsdata.track) <- dpars.data
      countsdata.track
    })
  }else{
    countsdata.tracks <- NULL
  }

  # Peak annotation tracks
  if(!missing(peaks) && (length(peaks) > 0)){
    dpars.peaks <- list(col.title = "black",
                        col.axis = "black",
                        col.border.title = "lightgray",
                        col.frame = "lightgray",
                        rotation.title = rotation.title,
                        frame = FALSE,
                        cex.axis = 0.2)

    peaks.tracks <- lapply(names(peaks), function(x){
      if(verbose){ cat("Adding", x, "track...\n") }
      peaks.gr <- peaks[[x]]
      seqlevelsStyle(peaks.gr) <- "UCSC"
      seqlevels(peaks.gr, pruning.mode = "coarse") <- paste0("chr",1:22)
      peaks.gr$score <- 1
      peaks.gr <- peaks.gr[countOverlaps(peaks.gr, region) > 0]
      peaks.track <- DataTrack(range = peaks.gr,
                               type = "h",
                               genome = genome,
                               name = x,
                               col = "navy", showAxis=F,
                               ylim = c(0,1))
      displayPars(peaks.track) <- dpars.peaks
      peaks.track
    })
  }else{
    peaks.tracks <- NULL
  }

  # HiC loops tracks
  if(!missing(HiC_loops) && (length(HiC_loops) > 0)){
    dpars.HiC <- list(col.interactions = "black",
                      col.anchors.fill = "blue",
                      col.anchors.line = "black",
                      interaction.dimension = "width",
                      interaction.measure = "score",
                      plot.trans = FALSE,
                      plot.outside = FALSE,
                      col.outside="lightblue",
                      anchor.height = 0.1,
                      rotation.title = rotation.title,
                      col.title = "black",
                      col.axis = "black",
                      col.border.title = "lightgray",
                      col.frame = "lightgray",
                      frame = FALSE,
                      cex.axis = 0.2)

    gene.annots$chr <- as.character(seqnames(gene.annots))
    gene.annots$tss <- start(resize(gene.annots, width = 1))

    HiC_loops.tracks <- lapply(names(HiC_loops), function(x){
      if(verbose){ cat("Adding", x, "track...\n") }
      HiC_loops.gr <- HiC_loops[[x]]
      if(filter_protein_coding_genes){
        HiC_loops.gr <- HiC_loops.gr[HiC_loops.gr$gene_name %in% gene.annots$gene_name]
      }

      HiC_promoters.gr <- GRanges(seqnames = HiC_loops.gr$promoter_chr,
                                  ranges = IRanges(start = HiC_loops.gr$promoter_start, end = HiC_loops.gr$promoter_end),
                                  score = HiC_loops.gr$score,
                                  gene = HiC_loops.gr$gene_name)
      HiC_enhancers.gr <- GRanges(seqnames = seqnames(HiC_loops.gr),
                                  ranges = IRanges(start = start(HiC_loops.gr), end = end(HiC_loops.gr)))
      HiC_loops.obj <- GenomicInteractions::GenomicInteractions(anchor1 = HiC_promoters.gr, anchor2 = HiC_enhancers.gr)
      HiC_loops.obj$counts <- round(HiC_loops.obj$anchor1.score)

      if(!is.null(filter_HiCloops_genes)){
        cat("Only shows", x, "links to", paste(filter_HiCloops_genes, collapse = ","), "\n")
        HiC_loops.obj <- HiC_loops.obj[which(HiC_loops.obj$anchor1.gene %in% filter_HiCloops_genes),]
      }

      if(!is.null(filter_HiCloops_snps)){
        cat("Only shows", x, "links to", paste(filter_HiCloops_snps, collapse = ","), "\n")
        highlighted.snps.gr <- finemapstats[finemapstats$snp %in% filter_HiCloops_snps]
        HiC_loops.obj <- subsetByOverlaps(HiC_loops.obj, highlighted.snps.gr)
      }

      HiC_loops.track <- InteractionTrack(HiC_loops.obj, name = x)
      displayPars(HiC_loops.track) <- dpars.HiC
      HiC_loops.track
    })
  }else{
    HiC_loops.tracks <- NULL
  }

  # gene track
  gene.track <- make_genetrack_obj(region, genetrack_db, txdb, gene.annots, genome, track_name = "Genes")

  dpars.genes <- list(col.title = "black",
                      col.axis = "black",
                      col.border.title = "lightgray",
                      col.frame = "lightgray",
                      rotation.title = rotation.title,
                      frame = FALSE,
                      cex.axis = 0.2)
  displayPars(gene.track) <- dpars.genes

  # # restrict to protein coding genes
  if(filter_protein_coding_genes){
    gene.track <- gene.track[symbol(gene.track) %in% gene.annots$gene_name]
  }

  # axis track
  axisTrack <- GenomeAxisTrack()

  # List all tracks
  list.of.tracks <- c(pval.track,
                      pip.track,
                      countsdata.tracks,
                      peaks.tracks,
                      HiC_loops.tracks,
                      gene.track,
                      axisTrack)
  track.sizes <- c(1,
                   0.6,
                   rep(0.4, length(countsdata.tracks)),
                   rep(0.4, length(peaks.tracks)),
                   rep(0.6, length(HiC_loops.tracks)),
                   0.5,
                   0.4)

  # Highlight SNPs
  if( !missing(highlight_snps) && (length(highlight_snps) > 0)){

    if("topSNP" %in% highlight_snps){
      topsnp_idx <- which.max(curr_finemapstats$pip)
      highlight_snps[which(highlight_snps == "topSNP")] <- curr_finemapstats$snp[topsnp_idx]
    }
    highlight_snps <- intersect(highlight_snps, curr_finemapstats$snp)
    highlight_pos <- curr_finemapstats$pos[match(highlight_snps, curr_finemapstats$snp)]
    cat("Highlight SNPs:", highlight_snps[which(highlight_snps %in% curr_finemapstats$snp)], "\n")

    if(length(highlight_pos) == 1){
      list.of.tracks <- HighlightTrack(trackList = list.of.tracks,
                                       start = c(highlight_pos-500), width = 1000,
                                       chromosome = as.character(seqnames(region)),
                                       col = highlight_colors)
    }else if(length(highlight_pos) > 1){
      cat("Highlight SNP positions:", highlight_pos, "\n")

      if(length(highlight_colors) > 1){
        highlight_colors <- highlight_colors[which(highlight_snps %in% curr_finemapstats$snp)]
      }
      list.of.tracks <- HighlightTrack(trackList = list.of.tracks,
                                       start = highlight_pos, width = 1,
                                       chromosome = as.character(seqnames(region)),
                                       col = highlight_colors)
    }
  }

  plotTracks(list.of.tracks,
             chromosome = as.character(seqnames(region)),
             transcriptAnnotation = "symbol",
             collapseTranscripts= "longest",
             from = start(region),
             to = end(region),
             sizes = track.sizes,
             just.group = genelabel_side)

}


# Generating gene track object for Gviz
make_genetrack_obj <- function(curr.locus.gr,
                               genetrack_db = c("txdb", "gene.annots","UCSC"),
                               txdb = NULL, gene.annots = NULL,
                               genome = "hg19", keytype = "ENSEMBL",
                               track_name = "Genes",
                               filter_protein_coding_genes = TRUE){

  genetrack_db <- match.arg(genetrack_db)

  if(genetrack_db == "txdb"){
    # use supplied txdb if available.
    cat("Making gene track object using gene annotations in txdb ...\n")
    grtrack <- GeneRegionTrack(range = txdb,
                               genome = genome,
                               chromosome = as.character(seqnames(curr.locus.gr)),
                               start = start(curr.locus.gr),
                               end = end(curr.locus.gr),
                               name = track_name)

    symbol(grtrack) <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                              keys=sub("\\.\\d.*$", "", gene(grtrack)),
                              keytype=keytype, column="SYMBOL")
    # symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), gene(grtrack), symbol(grtrack))
    symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), "", symbol(grtrack))

  }else if(genetrack_db == "gene.annots"){
    # use supplied gene.annots if available.
    cat("Making gene track object using gene annotations in gene.annots ...\n")
    grtrack <- GeneRegionTrack(range = gene.annots,
                               genome = genome,
                               chromosome = as.character(seqnames(curr.locus.gr)),
                               start = start(curr.locus.gr),
                               end = end(curr.locus.gr),
                               name = track_name,
                               symbol = gene.annots$gene_name)

  }else if(genetrack_db == "UCSC"){
    if(genome == "hg19"){
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      cat("Making gene track object using TxDb.Hsapiens.UCSC.hg19.knownGene ...\n")
    }else if(genome == "hg38"){
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      cat("Making gene track object using TxDb.Hsapiens.UCSC.hg38.knownGene ...\n")
    }else{
      stop("Genome only supports hg19 or hg38! Otherwise, please specify txdb or gene.annots!\n")
    }

    keytype <- "ENTREZID"

    grtrack <- GeneRegionTrack(range = txdb,
                               genome = genome,
                               chromosome = as.character(seqnames(curr.locus.gr)),
                               start = start(curr.locus.gr),
                               end = end(curr.locus.gr),
                               name = track_name)

    symbol(grtrack) <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                              keys=sub("\\.\\d+$", "", gene(grtrack)),
                              keytype=keytype, column="SYMBOL")
    symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), gene(grtrack), symbol(grtrack))
  }


  # restrict to protein coding genes
  if(filter_protein_coding_genes){
    grtrack <- grtrack[symbol(grtrack) %in% gene.annots$gene_name]
  }

  return(grtrack)
}

