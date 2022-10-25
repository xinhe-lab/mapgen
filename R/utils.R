
#' @title Process pcHiC data and save as a GRanges object
#'
#' @param pcHiC a data frame of pcHiC data, with columns named "Promoter" and
#' "Interacting_fragment". Interacting_fragment should contains
#' chr, start and end positions of the fragments interacting with promoters
#' e.g. "chr.start.end" or "chr:start-end".
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @return a GRanges object with processed pcHiC links, with genomic coordinates
#' of the interacting regions and gene names (promoters).
#' @export
process_pcHiC <- function(pcHiC){

  pcHiC <- pcHiC %>% dplyr::select(Promoter, Interacting_fragment)
  # separate genes connecting to the same fragment
  pcHiC <- pcHiC %>% tidyr::separate_rows(Promoter) %>% dplyr::rename(gene_name = Promoter)

  pcHiC <- pcHiC %>%
    tidyr::separate(Interacting_fragment, c('otherEnd_chr', 'otherEnd_start', 'otherEnd_end')) %>%
    dplyr::mutate(otherEnd_start = as.numeric(otherEnd_start), otherEnd_end = as.numeric(otherEnd_end))

  pcHiC.gr <- makeGRangesFromDataFrame(pcHiC,
                                       seqnames.field = 'otherEnd_chr',
                                       start.field = 'otherEnd_start',
                                       end.field = 'otherEnd_end',
                                       keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(pcHiC.gr) <- 'UCSC'

  return(pcHiC.gr)
}


#' @title Process ABC scores and save as a GRanges object
#'
#' @param ABC a data frame of ABC scores from Nasser et al. Nature 2021 paper
#' @param ABC.thresh Significance threshold of ABC scores
#' (default = 0.015, as in Nasser et al. Nature 2021 paper).
#' @param full.element Logical; if TRUE, use full length of ABC elements
#' extracted from the "name" column. Otherwise, use the original (narrow)
#' regions provided in the ABC scores data.
#' @param flank  Flanking regions around ABC elements (default = 0).
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @return a GRanges object with processed ABC scores, with genomic coordinates
#' of the interacting regions and gene names (promoters).
#' @export
process_ABC <- function(ABC, ABC.thresh = 0.015, full.element = FALSE, flank = 0){

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

  if(flank > 0){
    ABC$start <- ABC$start - flank
    ABC$end <- ABC$end + flank
  }

  ABC.gr <- makeGRangesFromDataFrame(ABC, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(ABC.gr) <- 'UCSC'

  return(ABC.gr)
}

#' @title Load ENCODE narrow peak data and convert to a GRanges object
#'
#' @param peak.file ENCODE narrow peak file
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @return a GRanges object for the peaks
#' @export
#'
process_narrowpeaks <- function(peak.file){
  if( !file.exists(peak.file) ){stop('narrowpeak file is not availble!')}
  peaks <- data.table::fread(peak.file)
  colnames(peaks) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')
  peaks.gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(peaks.gr) <- 'UCSC'
  return(peaks.gr)
}

# helper function to make a GRanges object
make_ranges <- function(seqname, start, end){
  return(GenomicRanges::GRanges(seqnames = seqname,
                                ranges = IRanges::IRanges(start = start, end = end)))
}


# get LD (r^2) between each SNP and the top SNP in sumstats
get_LD_bigSNP <- function(sumstats, bigSNP, topSNP = NULL){

  # only include SNPs in bigSNP markers
  sumstats <- sumstats[sumstats$snp %in% bigSNP$map$marker.ID, ]
  sumstats$bigSNP_idx <- match(sumstats$snp, bigSNP$map$marker.ID)

  if(missing(topSNP)){
    if( min(sumstats$pval) >=0 && max(sumstats$pval) <= 1 ){
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

