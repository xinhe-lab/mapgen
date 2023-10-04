
#' @title Process PC-HiC data and save as a GRanges object
#'
#' @param pcHiC A data frame of PC-HiC data, with columns named "Promoter" and
#' "Interacting_fragment". Interacting_fragment should contains
#' chr, start and end positions of the fragments interacting with promoters
#' e.g. "chr.start.end" or "chr:start-end".
#' @import GenomicRanges
#' @import tidyverse
#' @return A GRanges object with processed PC-HiC links, with genomic coordinates
#' of the interacting regions and gene names (promoters).
#' @export
process_pcHiC <- function(pcHiC){

  pcHiC <- pcHiC %>% dplyr::select(Promoter, Interacting_fragment)
  # separate genes connecting to the same fragment
  pcHiC <- pcHiC %>%
    tidyr::separate_rows(Promoter) %>%
    dplyr::rename(gene_name = Promoter)

  pcHiC <- pcHiC %>%
    tidyr::separate(Interacting_fragment, c('otherEnd_chr', 'otherEnd_start', 'otherEnd_end')) %>%
    dplyr::mutate(otherEnd_start = as.numeric(otherEnd_start), otherEnd_end = as.numeric(otherEnd_end))

  pcHiC.gr <- GenomicRanges::makeGRangesFromDataFrame(pcHiC,
                                       seqnames.field = 'otherEnd_chr',
                                       start.field = 'otherEnd_start',
                                       end.field = 'otherEnd_end',
                                       keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(pcHiC.gr) <- 'UCSC'

  return(pcHiC.gr)
}


#' @title Process ABC scores and save as a GRanges object
#'
#' @param ABC A data frame of ABC scores from Nasser et al. Nature 2021 paper
#' @param ABC.thresh Threshold of ABC scores
#' (default = 0.015, as in Nasser et al. Nature 2021 paper).
#' @param full.element Logical; if TRUE, use full length of ABC elements
#' extracted from the "name" column. Otherwise, use the original (narrow)
#' regions provided in the ABC scores data.
#' @param expand  Expand the ABC regions around ABC elements (default = 0).
#' @import GenomicRanges
#' @import tidyverse
#' @return a GRanges object with processed ABC scores, with genomic coordinates
#' of the interacting regions and gene names (promoters).
#' @export
process_ABC <- function(ABC, ABC.thresh = 0.015, full.element = FALSE, expand = 0){

  if(full.element){
    ABC <- ABC %>%
      tidyr::separate(name, c(NA, 'element_region'), sep = '\\|', remove = FALSE) %>%
      tidyr::separate(element_region, c(NA, 'element_location'), sep = '\\:') %>%
      tidyr::separate(element_location, c('element_start', 'element_end'), sep = '\\-') %>%
      dplyr::mutate(start = as.numeric(element_start), end = as.numeric(element_end))
  }

  ABC <- ABC %>%
    dplyr::rename(gene_name = TargetGene) %>%
    dplyr::filter(ABC.Score >= ABC.thresh)

  if(expand > 0){
    ABC$start <- ABC$start - expand
    ABC$end <- ABC$end + expand
  }

  ABC.gr <- GenomicRanges::makeGRangesFromDataFrame(ABC, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(ABC.gr) <- 'UCSC'

  return(ABC.gr)
}

# Load ENCODE narrow peak data and convert to a GRanges object
process_narrowpeaks <- function(peak.file){
  if( !file.exists(peak.file) ){stop('narrowpeak file is not availble!')}
  peaks <- data.table::fread(peak.file)
  colnames(peaks) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(peaks.gr) <- 'UCSC'
  return(peaks.gr)
}


# Annotations for causal SNPs
annotator_merged <- function(sumstats, annotations){

  snpRanges <- make_ranges(sumstats$chr, sumstats$pos, sumstats$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=sumstats$snp)
  sumstats['annots'] <- ''

  for(f in annotations){

    curr <- rtracklayer::import(f, format='bed')
    snpRangesIn <- IRanges::subsetByOverlaps(snpRanges, annot.gr)
    snpsIn <- unique(snpRangesIn$snp)

    if(length(snpsIn)>0){
      curr <- sumstats %>% dplyr::pull(annots)
      curr <- curr[sumstats$snp %in% snpsIn]
      delims <- rep(';', length(curr))
      delims[which(curr == '')] <- ''
      sumstats[sumstats$snp %in% snpsIn,"annots"] <-
        paste0(curr,delims,gsub(pattern = '.bed',replacement = '', x = basename(f)))
    }
  }
  return(sumstats)
}

# helper function to make a GRanges object
make_ranges <- function(seqname, start, end){
  return(GenomicRanges::GRanges(seqnames = seqname,
                                ranges = IRanges::IRanges(start = start, end = end)))
}

# get TSS positions from promoters or gene annotations.
get_tss <- function(gr, type = c('promoter', 'gene')){
  type <- match.arg(type)

  if(type == 'gene'){
    tss <- GenomicRanges::start(GenomicRanges::resize(gr, width = 1))
  }else if(type == 'promoter'){
    tss <- integer(length(gr))
    plus_strand <- which(as.character(GenomicRanges::strand(gr)) == '+')
    minus_strand <- which(as.character(GenomicRanges::strand(gr)) == '-')
    tss[plus_strand] <- GenomicRanges::end(gr)[plus_strand]
    tss[minus_strand] <- GenomicRanges::start(gr)[minus_strand]
  }

  return(tss)
}

# get LD (r^2) between each SNP and the top SNP in sumstats
get_LD_bigSNP <- function(sumstats, bigSNP, topSNP = NULL){

  # only include SNPs in bigSNP markers
  sumstats <- sumstats[sumstats$snp %in% bigSNP$map$marker.ID, ]
  sumstats$bigSNP_idx <- match(sumstats$snp, bigSNP$map$marker.ID)

  if(missing(topSNP)){
    if( max(sumstats$pval) <= 1 ){
      sumstats$pval <- -log10(sumstats$pval)
    }
    top_snp_idx <- sumstats$bigSNP_idx[which.max(sumstats$pval)]
  }else{
    top_snp_idx <- sumstats$bigSNP_idx[sumstats$snp == topSNP]
  }

  top_snp_genotype <- bigSNP$genotypes[,top_snp_idx]
  genotype.mat <- bigSNP$genotypes[,sumstats$bigSNP_idx]

  r2.vals <- as.vector(cor(top_snp_genotype, genotype.mat))^2
  sumstats$r2 <- round(r2.vals, 4)

  return(sumstats)
}


# Add LD information
add_LD_bigSNP <- function(sumstats, bigSNP,
                          r2.breaks = c(0, 0.1, 0.25, 0.75, 0.9, 1),
                          r2.labels = c("0-0.1","0.1-0.25","0.25-0.75","0.75-0.9","0.9-1")) {

  # only include SNPs in bigSNP markers
  sumstats <- sumstats[sumstats$snp %in% bigSNP$map$marker.ID, ]
  sumstats$bigSNP_idx <- match(sumstats$snp, bigSNP$map$marker.ID)

  if( max(sumstats$pval) <= 1 ){
    sumstats$pval <- -log10(sumstats$pval)
  }

  locus_list <- unique(sumstats$locus)

  sumstats.r2.df <- data.frame()
  for(locus in locus_list){
    curr_sumstats <- sumstats[sumstats$locus == locus, ]
    top_snp_idx <- curr_sumstats$bigSNP_idx[which.max(curr_sumstats$pval)]
    top_snp_genotype <- bigSNP$genotypes[,top_snp_idx]
    genotype.mat <- bigSNP$genotypes[,curr_sumstats$bigSNP_idx]

    r2.vals <- as.vector(cor(top_snp_genotype, genotype.mat))^2
    r2.brackets <- cut(r2.vals, breaks = r2.breaks, labels = r2.labels)
    curr_sumstats$r2 <- r2.brackets
    sumstats.r2.df <- rbind(sumstats.r2.df, curr_sumstats)
  }

  return(sumstats.r2.df)
}

# get gene region
get_gene_region <- function(gene.mapping.res, genes.of.interest, ext = 10000,
                            select.region = c("all", "locus")){

  select.region <- match.arg(select.region)

  high.conf.snp.df <- gene.mapping.res %>% dplyr::filter(pip > 0.2) %>%
    dplyr::group_by(snp) %>% dplyr::arrange(-gene_pip) %>% dplyr::slice(1)
  gene.gr <- gene.annots[match(high.conf.snp.df$gene_name, gene.annots$gene_name),]
  gene.gr$tss <- start(resize(gene.gr, width = 1))
  gene.gr <- gene.gr[,c("gene_name","tss")]
  high.conf.snp.df$tss <- gene.gr$tss

  gene.snp.tss <- high.conf.snp.df %>%
    dplyr::filter(gene_name %in% genes.of.interest) %>%
    dplyr::group_by(locus) %>%
    dplyr::arrange(-pip) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(distToTSS = pos-tss) %>%
    dplyr::select(gene_name, locus, chr, pos, tss, distToTSS)

  if(select.region == "all"){
    chr <- gene.snp.tss$chr[1]
    locus <- paste(gene.snp.tss$locus, collapse = ",")
    region_start <- min(c(gene.snp.tss$pos, gene.snp.tss$tss)) - ext
    region_end <- max(c(gene.snp.tss$pos, gene.snp.tss$tss)) + ext
    region <- GRanges(seqnames = chr,
                      IRanges(start = region_start, end = region_end),
                      locus = locus)
    seqlevelsStyle(region) <- "UCSC"
  }else if(select.region == "locus"){
    region <- lapply(gene.snp.tss$locus, function(l){
      locus <- l
      locus.snp.tss <- gene.snp.tss[gene.snp.tss$locus == locus, ]
      chr <- locus.snp.tss$chr
      if(distToTSS < 0){ #snp is upstream
        region_start <- locus.snp.tss$pos - ext
        region_end <- locus.snp.tss$tss + ext
      } else{ # snp is downstream
        region_start <- locus.snp.tss$tss - ext
        region_end <- locus.snp.tss$pos + ext
      }
      gene_locus_region.gr <- GRanges(seqnames = chr,
                                      IRanges(start = region_start, end = region_end),
                                      locus = locus)
      seqlevelsStyle(gene_locus_region.gr) <- "UCSC"
      gene_locus_region.gr
    })
  }
  return(region)
}
