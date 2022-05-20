
#' @title Prepare SuSiE summary statistics with TORUS prior probabilities
#' @description Adds TORUS SNP level prior probabilities to the summary statistics
#' @param sumstats a tibble or data frame containing raw summary statistics
#' @param torus_prior a tibble containing SNP level priors
#' (result from run_torus with \code{option=\dQuote{est-prior}})
#' @param torus_fdr a tibble containing the FDR of each region
#' (result from run_torus with \code{option=\dQuote{fdr}}).
#' Optional, if available, only keep the loci with Torus FDR < fdr_thresh.
#' @param fdr_thresh FDR cutoff (default: 0.1)
#' @return tibble of summary statistics updated with TORUS prior probabilities
#' @export
prepare_susie_data_with_torus_result <- function(sumstats, torus_prior, torus_fdr, fdr_thresh=0.1){

  if(!missing(torus_fdr)){
    # keep loci at fdr_thresh FDR (10% by default)
    chunks <- torus_fdr$region_id[torus_fdr$fdr < fdr_thresh]
    sumstats <- sumstats[sumstats$locus %in% chunks, ]
  }

  # Add Torus SNP-level prior probabilities
  sumstats <- dplyr::inner_join(sumstats, torus_prior, by='snp')

  return(sumstats)

}

#' @title Run fine-mapping
#' @description Run finemapping with SuSiE for all LD blocks
#' @param sumstats a tibble or data frame containing raw summary statistics; must have header!
#' @param bigSNP a bigsnpr object attached via bigsnpr::snp_attach()
#' @param priortype prior type: "torus" or "uniform".
#' @param L Number of causal signals
#' @return list of finemapping results; one per LD block
#' @export
run_finemapping <- function(sumstats, bigSNP, priortype = c('torus', 'uniform'), L = 1){

  priortype <- match.arg(priortype)

  if(priortype == 'torus'){
    useprior <- TRUE
    stopifnot('torus_prior' %in% colnames(sumstats))
  }else if(priortype == 'uniform'){
    useprior <- FALSE
  }

  chunks <- unique(sumstats$locus)
  susie_res <- list()
  for(i in seq_along(chunks)){
    z <- chunks[i]
    cat(sprintf('\rFinemapping locus %s...', z))
    susie.df <- sumstats[sumstats$locus == z, ]
    susie_res[[as.character(z)]] <- run_susie(susie.df, bigSNP, z, L, useprior)
    cat(sprintf('%.0f%% completed.', i/length(chunks)*100))
  }

  return(susie_res)

}


#' @title Run fine-mapping with SUSIE for one LD block
#' @param sumstats summary statistics
#' @param bigSNP bigSNP object
#' @param locus LD block index
#' @param L Number of causal signals
#' @param useprior Logical, if TRUE, use the \code{torus_prior} column
#' in \code{sumstats} as prior.
#' @return finemapping results for one LD block
#' @export
run_susie <- function(sumstats, bigSNP, locus, L=1, useprior){

  sub.sumstats <- sumstats[sumstats$locus == locus, ]
  if(nrow(sub.sumstats) > 1){
    X <- bigSNP$genotypes[ , sub.sumstats$bigSNP_index]
    X <- scale(X, center = T, scale = T)
    zhat <- sub.sumstats$zscore
    R <- cov2cor((crossprod(X) + tcrossprod(zhat))/nrow(X))
    if(useprior){
      res <- suppressWarnings(susieR::susie_rss(z = zhat,
                                                prior_weights = sub.sumstats$torus_prior,
                                                R = R,
                                                L = L,
                                                verbose = F))
    }
    else{
      res <- suppressWarnings(susieR::susie_rss(z = zhat,
                                                R = R,
                                                L = L,
                                                verbose = F))
    }
    return(res)
  }
}

#' @title merges SuSiE results with original summary statistics data frame
#' @description  merges SuSiE results with original summary statistics data frame
#' This function assumes L = 1. ONLY ONE CREDIBLE SET PER LOCUS!
#' @param susie_results data frame containing SuSiE finemapping result
#' @param sumstats data frame containing summary statistics
#'
#' @export
merge_susie_sumstats <- function(susie_results, sumstats){

  sumstats$susie_pip <- 0
  sumstats$CS <- 0
  loci <- names(susie_results)

  for(l in loci){
    n.snps <- length(susie_results[[l]]$pip)
    sumstats[sumstats$locus == as.numeric(l), 'susie_pip'] <- susie_results[[l]]$pip

    snps.in.cs <- rep(0, n.snps)
    if(!is.null(susie_results[[l]]$sets$cs)){
      snps.in.cs[unlist(susie_results[[l]]$sets$cs$L1)] <- 1
    }
    sumstats[sumstats$locus == as.numeric(l), 'CS'] <- snps.in.cs
  }

  return(sumstats)
}


#' @title Process fine-mapping summary statistics data
#'
#' @param finemapstats A data frame of fine-mapping summary statistics
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
#' @param maxL Maximum number of credible sets (default = 10).
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @return A GRanges object with cleaned and filtered fine-mapping summary statistics
#' @export
process_finemapping_sumstats <- function(finemapstats,
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
                                         maxL = 10,
                                         cols.to.keep = c('snp','chr','pos', 'pip', 'pval', 'zscore','cs', 'locus')){

  cat('Processing fine-mapping summary statistics ...\n')
  finemapstats <- finemapstats %>% dplyr::rename(snp = all_of(snp),
                                                 chr = all_of(chr),
                                                 pos = all_of(pos),
                                                 pip = all_of(pip))

  if( pval %in% colnames(finemapstats) ){
    finemapstats <- dplyr::rename(finemapstats, pval = all_of(pval))
  }else{
    finemapstats$pval <- NA
  }

  if( zscore %in% colnames(finemapstats) ){
    finemapstats <- dplyr::rename(finemapstats, zscore = all_of(zscore))
  }else{
    finemapstats$zscore <- NA
  }

  if( cs %in% colnames(finemapstats) ){
    finemapstats <- dplyr::rename(finemapstats, cs = all_of(cs))
  }else{
    finemapstats$cs <- NA
  }

  if( locus %in% colnames(finemapstats) ){
    finemapstats <- dplyr::rename(finemapstats, locus = all_of(locus))
  }else{
    finemapstats$locus <- NA
  }

  # Remove SNPs with multiple PIPs
  if(any(duplicated(paste(finemapstats$chr, finemapstats$pos)))){
    cat('Remove SNPs with multiple PIPs...\n')
    finemapstats <- finemapstats %>% dplyr::arrange(desc(pip)) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)
  }

  finemapstats.gr <- makeGRangesFromDataFrame(finemapstats, start.field = 'pos', end.field = 'pos', keep.extra.columns = TRUE)
  finemapstats.gr$chr <- finemapstats$chr
  finemapstats.gr$pos <- finemapstats$pos
  mcols(finemapstats.gr) <- mcols(finemapstats.gr)[,cols.to.keep]
  GenomeInfoDb::seqlevelsStyle(finemapstats.gr) <- 'UCSC'

  if( pip.thresh > 0 ) {
    cat('Filter SNPs with PIP threshold of', pip.thresh, '\n')
    finemapstats.gr <- finemapstats.gr[finemapstats.gr$pip > pip.thresh, ]
  }

  if( filterCS ) {
    cat('Filter SNPs in credible sets \n')
    finemapstats.gr <- finemapstats.gr[finemapstats.gr$cs >= 1 & finemapstats.gr$cs <= maxL, ]
  }

  return(finemapstats.gr)

}
