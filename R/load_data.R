
#' @title Load finemapping summary statistics data
#'
#' @param finemapping.file Filename of finemapping summary statistics data
#' @param SNP Name of the SNP ID (rsID) column in the finemapping summary statistics data
#' @param CHR Name of the CHR column in the finemapping summary statistics data
#' @param POS Name of the POS column in the finemapping summary statistics data
#' @param PVAL Name of the P-value column in the finemapping data frame
#' @param CS Name of the CS column in the finemapping data frame
#' @import tidyverse
#' @export
#'
load_finemapping_sumstat <- function(finemapping.file,
                                     SNP = 'SNP',
                                     CHR = 'CHR',
                                     POS = 'BP',
                                     PVAL = 'pval',
                                     CS = 'cs'){

  cat('Load fine-mapping summary statistics data ... \n')
  finemap <- data.table::fread(finemapping.file)

  finemap <- finemap %>% dplyr::rename(snp = all_of(SNP), chr = all_of(CHR), pos = all_of(POS))

  if( PVAL %in% colnames(finemap) ){
    finemap <- dplyr::rename(finemap, pval = all_of(PVAL))
  }else{
    finemap$pval <- NA
  }

  if( CS %in% colnames(finemap) ){
    finemap <- dplyr::rename(finemap, cs = all_of(CS))
  }else{
    finemap$cs <- NA
  }

  return(finemap)

}


#' @title Process fine-mapping summary statistics data
#'
#' @param finemap A data frame of fine-mapping summary statistics
#' @param pip.thresh PIP threshold (default = 1e-5).
#' @param filterCS If TRUE, limiting to SNPs within credible sets.
#' @param maxCS Maximum number of credible sets (default = 10).
#' @import tidyverse
#' @export
process_finemapping_sumstat <- function(finemap, pip.thresh = 1e-5, filterCS = FALSE, maxCS = 10){
  finemap <- finemap %>% dplyr::select(snp, pip, chr, pos, pval, cs)

  # Remove duplicated SNPs with multiple PIPs
  finemap <- finemap %>% dplyr::arrange(desc(pip)) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)

  finemap.gr <- GenomicRanges::makeGRangesFromDataFrame(finemap, start.field = 'pos', end.field = 'pos', keep.extra.columns = TRUE)
  finemap.gr$chr <- finemap$chr
  finemap.gr$pos <- finemap$pos
  GenomeInfoDb::seqlevelsStyle(finemap.gr) <- 'UCSC'

  if( pip.thresh > 0 ) {
    cat('Filter SNPs with pip >', pip.thresh, '\n')
    finemap.gr <- finemap.gr[finemap.gr$pip > pip.thresh, ]
  }

  if( filterCS ) {
    cat('Filter SNPs in credible sets \n')
    finemap.gr <- finemap.gr[finemap.gr$cs >= 1 & finemap.gr$cs <= maxCS, ]
  }

  return(finemap.gr)

}

#' @title Load genomic annotations
#'
#' @param genomic_annots_file Genomic annotation rds file
#'
#' @export
load_genomic_annots <- function(genomic.annots.file){

  if( !file.exists(genomic.annots.file) ){
    stop('Genomic annotation file is not availble!')
  }
  cat('Load genomic annotations ...\n')

  genomic.annots <- readRDS(genomic.annots.file)

  genomic.annots$promoters$tss <- NA
  plus_strand <- which(strand(genomic.annots$promoters) == '+')
  minus_strand <- which(strand(genomic.annots$promoters) == '-')
  genomic.annots$promoters$tss[plus_strand] <- end(genomic.annots$promoters[plus_strand])
  genomic.annots$promoters$tss[minus_strand] <- start(genomic.annots$promoters[minus_strand])

  return(genomic.annots)
}

#' @title Load pcHiC data and save as a GRanges object
#'
#' @param pcHiC.file pcHiC file downloaded from Jung et al.
#' (https://pubmed.ncbi.nlm.nih.gov/31501517/)
#' @import tidyverse
#' @export
load_pcHiC <- function(pcHiC.file){

  if( !file.exists(pcHiC.file) ){stop('pcHiC file is not availble!')}
  cat('Load pcHiC data...\n')

  pcHiC <- data.table::fread(pcHiC.file)

  pcHiC <- pcHiC %>% dplyr::select(Promoter, Interacting_fragment)
  # separate genes connecting to the same fragment
  pcHiC <- pcHiC %>% tidyr::separate_rows(Promoter) %>% dplyr::rename(gene_name = Promoter)

  pcHiC <- pcHiC %>%
    tidyr::separate(Interacting_fragment, c('otherEnd_chr', 'otherEnd_start', 'otherEnd_end'), '\\.') %>%
    dplyr::mutate(otherEnd_start = as.numeric(otherEnd_start), otherEnd_end = as.numeric(otherEnd_end))

  pcHiC.gr <- GenomicRanges::makeGRangesFromDataFrame(pcHiC,
                                       seqnames.field = 'otherEnd_chr',
                                       start.field = 'otherEnd_start',
                                       end.field = 'otherEnd_end',
                                       keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(pcHiC.gr) <- 'UCSC'

  return(pcHiC.gr)
}


#' @title Load ABC scores and save as a GRanges object
#'
#' @param ABC.file ABC file downloaded from Nasser et al.
#' (https://www.engreitzlab.org/resources).
#' @param ABC.thresh Significance threshold of ABC scores (default = 0).
#' @param full.element Logitic; if TRUE, use full length of ABC elements.
#' Otherwise, use the narrow regions provided in the ABC data.
#' @param flank  Flanking regions around ABC elements (default = 0).
#' @import tidyverse
#' @export
load_ABC <- function(ABC.file, ABC.thresh = 0, full.element = TRUE, flank = 0){

  if( !file.exists(ABC.file) ){stop('ABC file is not availble!')}

  cat('Load ABC data...\n')

  ABC <- data.table::fread(ABC.file)

  if(full.element){
    ABC <- ABC %>%
      tidyr::separate(name, c(NA, 'element_region'), sep = '\\|', remove = FALSE) %>%
      tidyr::separate(element_region, c(NA, 'element_location'), sep = '\\:') %>%
      tidyr::separate(element_location, c('element_start', 'element_end'), sep = '\\-') %>%
      tidyr::mutate(start = as.numeric(element_start), end = as.numeric(element_end))
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
#' @export
#'
load_narrowpeaks <- function(peak.file){
  if( !file.exists(ABC.file) ){stop('ABC file is not availble!')}
  peaks <- data.table::fread(peak.file)
  colnames(peaks) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(peaks.gr) <- 'UCSC'
  return(peaks.gr)
}

#' @title Filtering coaccessibility data
#'
#' @param coaccess data frame of coaccessibility
#' @param cor.thresh Threshold for coaccessibility
#' @param dist.thresh Threshold for distance
#'
#' @export
#'
filter_coaccess <- function(coaccess, cor.thresh = 0.7, dist.thresh = 1e6){
  coaccess$dist <- abs(round((start(coaccess)+end(coaccess))/2) - round((coaccess$promoter_start + coaccess$promoter_end)/2))
  coaccess <- coaccess[coaccess$correlation > cor.thresh & coaccess$dist < dist.thresh,]
  return(coaccess)
}
