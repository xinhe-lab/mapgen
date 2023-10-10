
#' @title Prepare summary statistics with TORUS SNP-level priors as
#' input data for SuSiE
#' @description Check for required columns in summary statistics and
#' adds TORUS SNP-level priors to summary statistics to be used as input
#' data for \code{susie_rss}.
#' @param sumstats A data frame of summary statistics
#' @param torus_prior A data frame with SNP level priors
#' (result from \code{run_torus()} with \code{option=\dQuote{est-prior}})
#' @param torus_fdr A data frame containing the FDR of each region
#' (result from \code{run_torus()} with \code{option=\dQuote{fdr}}).
#' Optional, if available, only keep the loci with FDR < `fdr.thresh`.
#' @param fdr.thresh FDR cutoff (default: 0.1)
#' @return A data frame of summary statistics with SNP-level priors
#' to be used as input data for \code{susie_rss}.
#' @export
prepare_susie_data_with_torus_result <- function(sumstats,
                                                 torus_prior,
                                                 torus_fdr,
                                                 fdr.thresh = 0.1){

  # Check for required columns in sumstats
  required.cols <- c('chr','pos','snp','pval','locus','bigSNP_index')

  if(!all(required.cols %in% colnames(sumstats))){
    stop(sprintf('Column \"%s\" cannot be found in the summary statistics!',
                 required.cols[which(!required.cols %in% colnames(sumstats))]))
  }

  if(!'zscore' %in% colnames(sumstats)){
    cat("'zscore' not found in sumstats. Computing z-scores using beta and se ...\n")
    sumstats$zscore <- sumstats$beta / sumstats$se
  }

  if(!missing(torus_fdr)){
    # Select loci by FDR from TORUS
    cat('Select loci by TORUS FDR <', fdr.thresh, '\n')
    selected.loci <- unique(torus_fdr$region_id[torus_fdr$fdr < fdr.thresh])
    if(length(selected.loci) == 0){
      stop('No loci selected!')
    }else{
      sumstats <- sumstats[sumstats$locus %in% selected.loci, ]
    }
  }

  # Add SNP-level priors
  sumstats <- dplyr::inner_join(sumstats, torus_prior, by = 'snp')

  return(sumstats)
}

#' @title Run fine-mapping using summary statistics
#' @description Run fine-mapping with SuSiE using summary statistics
#' for all LD blocks with prior probabilities computed by TORUS
#' @param sumstats A data frame of summary statistics
#' @param bigSNP a \code{bigsnpr} object attached via \code{bigsnpr::snp_attach()}
#' containing the reference genotype panel.
#' @param LD_matrices A list of LD matrices (R, correlation matrices) for
#' all the LD blocks in `sumstats`, names of the list should
#' correspond to the 'locus' column in `sumstats`.
#' @param n The sample size (optional, but strongly recommended.)
#' @param priortype prior type:
#' 'torus' (use the 'torus_prior' in `sumstats`)
#' or 'uniform' (uniform prior).
#' @param L Number of causal signals.
#' @param estimate_residual_variance The default is FALSE,
#' the residual variance is fixed to 1 or variance of y.
#' If the in-sample LD matrix is provided,
#' we recommend setting estimate_residual_variance = TRUE.
#' @param verbose If TRUE, print progress,
#' and a summary of \code{susie_rss}.
#' @return A list of SuSiE results; one per LD block.
#' @export
run_finemapping <- function(sumstats,
                            bigSNP,
                            LD_matrices,
                            n,
                            priortype = c('torus', 'uniform'),
                            L = 1,
                            estimate_residual_variance = FALSE,
                            verbose = FALSE,
                            ...){

  priortype <- match.arg(priortype)

  if(priortype == 'torus'){
    useprior <- TRUE
    stopifnot('torus_prior' %in% colnames(sumstats))
  }else if(priortype == 'uniform'){
    useprior <- FALSE
  }

  finemap.locus.list <- unique(sumstats$locus)

  if(!missing(LD_matrices)){
    if(!setequal(names(LD_matrices), finemap.locus.list)){
      stop("Names in LD matrices do not match with the list of loci in sumstats!")
    }
  }

  susie.res <- list()
  for(locus in finemap.locus.list){
    cat(sprintf('Finemapping locus %s...\n', locus))
    sumstats_locus <- sumstats[sumstats$locus == locus, ]
    if(!missing(LD_matrices)){
      # use R from LD_matrices
      R <- LD_matrices[[as.character(locus)]]
      if(verbose){ cat("Using R from LD_matrices...\n")}
      susie.res[[as.character(locus)]] <- run_susie_rss(sumstats_locus,
                                                        R=R,
                                                        n=n,
                                                        L=L,
                                                        useprior=useprior,
                                                        estimate_residual_variance=estimate_residual_variance,
                                                        verbose=verbose,
                                                        ...)
    }else{
      # compute R using bigSNP reference genotype panel
      if(missing(bigSNP)){
        stop("Please provide LD matrix or bigSNP object!")
      }
      if(verbose){ cat("Computing R using bigSNP genotype matrix...\n") }
      susie.res[[as.character(locus)]] <- run_susie_rss(sumstats_locus,
                                                        bigSNP=bigSNP,
                                                        n=n,
                                                        L=L,
                                                        useprior=useprior,
                                                        estimate_residual_variance=estimate_residual_variance,
                                                        verbose=verbose,
                                                        ...)
    }

    if(verbose){
      cat(sprintf('%.0f%% completed.\n', length(susie.res)/length(finemap.locus.list)*100))
    }
  }

  return(susie.res)

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
    sumstats[sumstats$locus == as.numeric(l), 'cs'] <- snps.in.cs
  }

  return(sumstats)
}


#' @title Process fine-mapping summary statistics data
#'
#' @param finemapstats A data frame of fine-mapping summary statistics
#' @param snp Name of the SNP ID (rsID) column in the fine-mapping summary statistics
#' @param chr Name of the chr column in the fine-mapping summary statistics
#' @param pos Name of the position column in the fine-mapping summary statistics
#' @param pip Name of the PIP column in the fine-mapping summary statistics
#' @param pval Name of the P-value column in the fine-mapping summary statistics
#' @param zscore Name of the z-score column in the fine-mapping summary statistics
#' @param cs Name of the credible set (CS) column in the fine-mapping summary statistics
#' @param locus Name of the locus column in the fine-mapping summary statistics
#' @param pip.thresh Select SNPs by PIP threshold (default = 0, no filtering).
#' @param filterCS If TRUE, limiting to SNPs within credible sets.
#' @param maxL Maximum number of credible sets (default = 10).
#' If filterCS is TRUE, it will only keep SNPs with credible set (CS) number >= 1
#' and <= `maxL`.
#' @import GenomicRanges
#' @import tidyverse
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
                                         pip.thresh = 0,
                                         filterCS = FALSE,
                                         maxL = 10){

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
    finemapstats <- finemapstats %>% dplyr::arrange(desc(pip)) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)
  }

  finemapstats.gr <- GenomicRanges::makeGRangesFromDataFrame(finemapstats,
                                                             start.field = 'pos', end.field = 'pos',
                                                             keep.extra.columns = TRUE)

  finemapstats.gr <- plyranges::mutate(finemapstats.gr,
                                       chr=finemapstats$chr, pos = finemapstats$pos)

  GenomeInfoDb::seqlevelsStyle(finemapstats.gr) <- 'UCSC'

  if( pip.thresh > 0 ) {
    cat('Select SNPs with PIP >', pip.thresh, '\n')
    finemapstats.gr <- finemapstats.gr[finemapstats.gr$pip > pip.thresh, ]
  }

  if( filterCS ) {
    cat('Select SNPs within credible sets. \n')
    finemapstats.gr <- finemapstats.gr[finemapstats.gr$cs >= 1 & finemapstats.gr$cs <= maxL, ]
  }

  return(finemapstats.gr)

}


# Run fine-mapping with SuSiE using summary statistics (susie_rss)
run_susie_rss <- function(sumstats,
                          bigSNP,
                          R,
                          n,
                          L=1,
                          useprior=FALSE,
                          estimate_residual_variance=FALSE,
                          verbose=FALSE,
                          ...){

  if(is.null(sumstats$zscore)){
    cat("'zscore' not found in sumstats. Computing z-scores using beta and se ...\n")
    sumstats$zscore <- sumstats$beta / sumstats$se
  }

  z <- sumstats$zscore

  if(missing(R)){
    if(missing(bigSNP)){
      stop("Please provide LD matrix or bigSNP object!")
    }
    # compute R using reference panel in bigSNP object
    X <- bigSNP$genotypes[, sumstats$bigSNP_index]
    X <- scale(X, center = TRUE, scale = TRUE)
    R <- cor(X)
  }

  if(useprior){
    prior_weights <- sumstats$torus_prior
  }else{
    prior_weights <- NULL
  }

  cat('Run susie_rss...\n')

  if(verbose){
    cat(sprintf('n=%d, L=%d, useprior=%s, estimate_residual_variance=%s...\n',
                n, L, useprior, estimate_residual_variance))
  }

  res <- susieR::susie_rss(z = z,
                           R = R,
                           n = n,
                           L = L,
                           prior_weights = prior_weights,
                           estimate_residual_variance = estimate_residual_variance,
                           verbose = verbose,
                           ...)
  return(res)

}
