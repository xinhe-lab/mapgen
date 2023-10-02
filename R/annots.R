
#' @title Make genomic annotations from a GTF file
#'
#' @param gtf.file GTF file name
#' @param save If TRUE, save genomic annotations as a RDS file
#' @param outname Output file name
#' @import GenomicRanges
#' @import tidyverse
#' @export
make_genomic_annots <- function(gtf.file, save = FALSE, outname = NULL) {

  my.gtf <- rtracklayer::import(con = gtf.file, format = 'gtf')
  seqlevels(my.gtf, pruning.mode = 'coarse') <- paste0('chr',1:22)

  my.gtf.protein <- my.gtf[which(my.gtf$gene_type=='protein_coding'),]

  my.genes <- my.gtf.protein[my.gtf.protein$type=='gene',]

  canonical.transcripts <- my.gtf.protein %>%
    as_tibble() %>%
    dplyr::filter(type == 'transcript') %>%
    group_by(gene_id) %>%
    mutate(transLength = abs(end - start)) %>%
    arrange(-transLength) %>%
    dplyr::slice(1)

  canoncial.transcripts.str <- canonical.transcripts %>% .$transcript_id

  # keep only canonical transcripts
  my.gtf.protein.canonical <- my.gtf.protein[my.gtf.protein$transcript_id %in% canoncial.transcripts.str,]
  my.exons <- my.gtf.protein.canonical[my.gtf.protein.canonical$type=='exon',]
  my.UTR <- my.gtf.protein.canonical[my.gtf.protein.canonical$type=='UTR',]

  # get intron coordinates
  my.introns <- substract_exons(my.genes, my.exons)

  # get promoter coordinates
  my.genes.plus <- my.genes[strand(my.genes)=='+',]
  my.genes.neg <- my.genes[strand(my.genes)=='-',]
  my.promoters.plus <- GRanges(seqnames(my.genes.plus),
                               ranges = IRanges(start = start(my.genes.plus) - 2000,
                                                end = start(my.genes.plus)),
                               strand = strand(my.genes.plus))
  my.promoters.neg <- GRanges(seqnames(my.genes.neg),
                              ranges = IRanges(start = end(my.genes.neg),
                                               end = end(my.genes.neg) + 2000),
                              strand = strand(my.genes.neg))
  my.promoters <- append(my.promoters.plus, my.promoters.neg)
  my.promoters$gene_name <- c(my.genes.plus$gene_name, my.genes.neg$gene_name)

  exon.chr <- c(as.character(seqnames(my.exons)), as.character(seqnames(my.exons)))
  exon.pos <- c(start(my.exons), end(my.exons))
  exon.gene.name <- c(my.exons$gene_name, my.exons$gene_name)
  my.splice.junc <- GRanges(seqnames = exon.chr,
                            ranges = IRanges(start = exon.pos-100, end = exon.pos+100),
                            gene_name = exon.gene.name)

  annots <- list(
    genes = my.genes,
    exons = my.exons,
    introns = my.introns,
    UTRs = my.UTR,
    promoters = my.promoters,
    splice_junctions = my.splice.junc
  )

  if(save){
    if(!dir.exists(dirname(outname)))
      dir.create(dirname(outname), showWarnings = FALSE, recursive = TRUE)

    saveRDS(annots, outname)
  }

  return(annots)
}

# Subtracts all exon coordinates from all gene coordinates to get all intron coordinates
substract_exons <- function( genes.gr, exons.gr ) {
  # adapted from bosberg on Biostars. (https://www.biostars.org/p/489350/)
  # Subtract GRange object gr2 from gr1, but unlike setdiff, preserve individual ranges in gr1
  genes.df <- data.frame( seqnames=seqnames(genes.gr),
                          start=start(genes.gr)-1, end=end(genes.gr), # convert to 0-based start for bedtools
                          strand=strand(genes.gr),
                          gene_name = genes.gr$gene_name )
  exons.df <- data.frame( seqnames=seqnames(exons.gr),
                          start=start(exons.gr)-1,
                          end=end(exons.gr),
                          strand=strand(exons.gr) )
  result <- bedtoolsr::bt.subtract(genes.df, exons.df)

  if ( length(result)==0 ){
    # subtraction has left nothing remaining. Return empty GRanges obj.
    return( GRanges() )
  } else {
    colnames(result) <- colnames(genes.df)
    result$start <- result$start+1 # reset to 1-based notation consistent with GRanges
    return( GRanges(result) )
  }
}

