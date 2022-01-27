
#' @title Process fine mapping summary statistics data
#'
#' @param finemap A data frame of fine-mapping summary statistics
#' @param snp Name of the SNP ID (rsID) column in the summary statistics data
#' @param chr Name of the chr column in the summary statistics data frame
#' @param pos Name of the position column in the summary statistics data frame
#' @param pip Name of the PIP column in the summary statistics data frame
#' @param pval Name of the P-value column in the summary statistics data frame
#' @param zscore Name of the z-score column in the summary statistics data frame
#' @param cs Name of the CS column in the summary statistics data frame
#' @param locus Name of the locus column in the summary statistics data frame
#' @param cols.to.keep columns to keep in the returned data frame
#' @param pip.thresh PIP threshold (default = 1e-5).
#' @param filterCS If TRUE, limiting to SNPs within credible sets.
#' @param maxCS Maximum number of credible sets (default = 10).
#' @import tidyverse
#' @return A GRanges object with cleaned and filtered fine-mapping summary statistics
#' @export
process_finemapping_sumstat <- function(finemap,
                                        snp = 'snp',
                                        chr = 'chr',
                                        pos = 'pos',
                                        pip = 'pip',
                                        pval = 'pval',
                                        zscore = 'zscore',
                                        cs = 'cs',
                                        locus = 'locus',
                                        pip.thresh = 1e-5,
                                        filterCS = FALSE,
                                        maxCS = 10,
                                        cols.to.keep = c('snp','chr','pos', 'pip', 'pval', 'zscore','cs', 'locus')){

  cat('Process fine-mapping summary statistics ...\n')
  finemap <- finemap %>% dplyr::rename(snp = all_of(snp),
                                       chr = all_of(chr),
                                       pos = all_of(pos),
                                       pip = all_of(pip))

  if( pval %in% colnames(finemap) ){
    finemap <- dplyr::rename(finemap, pval = all_of(pval))
  }else{
    finemap$pval <- NA
  }

  if( zscore %in% colnames(finemap) ){
    finemap <- dplyr::rename(finemap, zscore = all_of(zscore))
  }else{
    finemap$zscore <- NA
  }

  if( cs %in% colnames(finemap) ){
    finemap <- dplyr::rename(finemap, cs = all_of(cs))
  }else{
    finemap$cs <- NA
  }

  if( locus %in% colnames(finemap) ){
    finemap <- dplyr::rename(finemap, locus = all_of(locus))
  }else{
    finemap$locus <- NA
  }

  # Remove SNPs with multiple PIPs
  if(any(duplicated(paste(finemap$chr, finemap$pos)))){
    cat('Remove SNPs with multiple PIPs...\n')
    finemap <- finemap %>% dplyr::arrange(desc(pip)) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)
  }

  finemap.gr <- GenomicRanges::makeGRangesFromDataFrame(finemap, start.field = 'pos', end.field = 'pos', keep.extra.columns = TRUE)
  finemap.gr$chr <- finemap$chr
  finemap.gr$pos <- finemap$pos
  mcols(finemap.gr) <- mcols(finemap.gr)[,cols.to.keep]
  GenomeInfoDb::seqlevelsStyle(finemap.gr) <- 'UCSC'

  if( pip.thresh > 0 ) {
    cat('Filter SNPs with PIP threshold of', pip.thresh, '\n')
    finemap.gr <- finemap.gr[finemap.gr$pip > pip.thresh, ]
  }

  if( filterCS ) {
    cat('Filter SNPs in credible sets \n')
    finemap.gr <- finemap.gr[finemap.gr$cs >= 1 & finemap.gr$cs <= maxCS, ]
  }

  return(finemap.gr)

}

#' @title Process pcHiC data and save as a GRanges object
#'
#' @param pcHiC a data frame of pcHiC data, with columns named "Promoter" and
#' "Interacting_fragment". Interacting_fragment should contains
#' chr, start and end positions of the fragments interacting with promoters
#' e.g. "chr.start.end" or "chr:start-end".
#' @import tidyverse
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
#' @param ABC a data frame of ABC scores from Nasser et al. Nature 2021 paper
#' @param ABC.thresh Significance threshold of ABC scores
#' (default = 0.015, as in Nasser et al. Nature 2021 paper).
#' @param full.element Logical; if TRUE, use full length of ABC elements
#' extracted from the "name" column. Otherwise, use the original (narrow)
#' regions provided in the ABC scores data.
#' @param flank  Flanking regions around ABC elements (default = 0).
#' @import tidyverse
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

  ABC.gr <- GenomicRanges::makeGRangesFromDataFrame(ABC, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(ABC.gr) <- 'UCSC'

  return(ABC.gr)
}

#' @title Load ENCODE narrow peak data and convert to a GRanges object
#'
#' @param peak.file ENCODE narrow peak file
#' @return a GRanges object for the peaks
#' @export
#'
process_narrowpeaks <- function(peak.file){
  if( !file.exists(peak.file) ){stop('narrowpeak file is not availble!')}
  peaks <- data.table::fread(peak.file)
  colnames(peaks) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(peaks.gr) <- 'UCSC'
  return(peaks.gr)
}

