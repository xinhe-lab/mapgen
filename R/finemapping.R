
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
#' @param filter filter SNPs from sumstats:
#'   'none': do not filter SNPs,
#'   'pval': filter SNPs with GWAS pval < \code{pval.thresh},
#'   'torus_fdr': filter SNPs with TORUS FDR < \code{torus_fdr}.
#' @param fdr.thresh FDR cutoff (default: 0.1)
#' @param pval.thresh GWAS pval cutoff (default: 5e-8)
#' @return A data frame of summary statistics with SNP-level priors
#' to be used as input data for finemapping with SuSiE.
#' @export
prepare_susie_data_with_torus_result <- function(sumstats,
                                                 torus_prior,
                                                 torus_fdr = NULL,
                                                 filter = c("none", "pval", "torus_fdr"),
                                                 fdr.thresh = 0.1,
                                                 pval.thresh = 5e-8){

  filter <- match.arg(filter)

  # Check for required columns in sumstats
  required.cols <- c('chr','pos','snp','pval','locus','bigSNP_index')

  if(!all(required.cols %in% colnames(sumstats))){
    stop(sprintf('Column \"%s\" cannot be found in the summary statistics!',
                 required.cols[which(!required.cols %in% colnames(sumstats))]))
  }

  if(!'zscore' %in% colnames(sumstats)){
    cat("'zscore' not found in sumstats. Computing z-scores ...\n")
    sumstats$zscore <- sumstats$beta / sumstats$se
  }

  # add torus priors
  sumstats <- dplyr::inner_join(sumstats, torus_prior, by = 'snp')

  if (filter == 'pval') {
    cat('Select loci by pval <', pval.thresh, '\n')
    if(max(sumstats$pval) > 1)
      pval <- 10^(-sumstats$pval)
    selected.loci <- unique(sumstats$locus[pval < pval.thresh])

    if(length(selected.loci) == 0)
      stop('No loci after filtering!')

    sumstats <- sumstats[sumstats$locus %in% selected.loci, , drop = FALSE]
  } else if (filter == 'torus_fdr') {
    # Select loci by FDR from TORUS
    cat('Select loci by TORUS FDR <', fdr.thresh, '\n')
    selected.loci <- unique(torus_fdr$region_id[torus_fdr$fdr < fdr.thresh])

    if(length(selected.loci) == 0)
      stop('No loci after filtering!')

    sumstats <- sumstats[sumstats$locus %in% selected.loci, , drop = FALSE]
  }

  return(sumstats)
}


#' @title Run fine-mapping for one region with SuSiE using summary statistics (susie_rss)
#' @description Run fine-mapping with SuSiE using summary statistics for one LD block
#' @param sumstats A data frame of summary statistics
#' @param R LD (correlation) matrix
#' @param snp_info Variant information for the LD matrix.
#' @param bigSNP a \code{bigsnpr} object attached via \code{bigsnpr::snp_attach()}
#' containing the reference genotype panel.
#' @param n GWAS sample size (optional, but strongly recommended.)
#' @param prior_weights A vector of prior probability for each SNP
#' @param L Number of causal signals.
#' @param estimate_residual_variance The default is FALSE,
#' the residual variance is fixed to 1 or variance of y.
#' If the in-sample LD matrix is provided,
#' we recommend setting estimate_residual_variance = TRUE.
#' @param verbose If TRUE, print progress
#' @param save If TRUE, save SuSiE result and LD (R) matrix.
#' @return a list of SuSiE results, and z-scores and LD (R) matrix used in SuSiE
#' @export
susie_finemap_region <- function(sumstats,
                                 R=NULL,
                                 snp_info=NULL,
                                 bigSNP=NULL,
                                 n=NULL,
                                 L=1,
                                 prior_weights=NULL,
                                 estimate_residual_variance=FALSE,
                                 verbose=FALSE,
                                 save = FALSE,
                                 ...){

  if (!is.null(R)) {
    # match the SNPs between GWAS sumstats with LD reference
    res <- match_gwas_LDREF(sumstats, R, snp_info)
    sumstats <- res$sumstats
    R <- res$R
    rm(res)
  } else if (!is.null(bigSNP)) {
    # compute LD using reference panel in bigSNP object
    if(is.null(sumstats$bigSNP_index))
      stop('bigSNP_index is missing in sumstats')
    cat("Compute LD matrix using bigSNP...\n")
    X <- bigSNP$genotypes[, sumstats$bigSNP_index]
    X <- scale(X, center = TRUE, scale = TRUE)
    R <- cor(X)
  } else if (L == 1) {
    # R does not matter for susie when L = 1
    R <- diag(nrow(sumstats))
  } else {
    stop("Please provide R (LD) matrix or bigSNP object when L > 1")
  }

  if (is.null(sumstats$zscore)) {
    cat("Compute z-scores\n")
    sumstats$zscore <- sumstats$beta / sumstats$se
  }

  z <- sumstats$zscore

  cat('Run susie_rss...\n')
  susie.res <- susieR::susie_rss(z = z,
                                 R = R,
                                 n = n,
                                 L = L,
                                 prior_weights = prior_weights,
                                 estimate_residual_variance = estimate_residual_variance,
                                 verbose = verbose,
                                 ...)

  susie.res$id <- sumstats$snp

  if (save) {
    susie.res$R <- R
    susie.res$z <- z
  }

  return(susie.res)

}

#' @title Run functional fine-mapping using summary statistics
#' @description Run fine-mapping with SuSiE using summary statistics
#' for all LD blocks with prior probabilities
#' @param sumstats A data frame of summary statistics
#' @param bigSNP a \code{bigsnpr} object attached via \code{bigsnpr::snp_attach()}
#' containing the reference genotype panel.
#' @param region_info A data frame of region information,
#' paths of LD matrices (R, correlation matrices),
#' and paths of variant information files corresponding to the LD matrices.
#' @param n The sample size (optional, but strongly recommended.)
#' @param priortype prior type:
#'   'torus' (use the 'torus_prior' in `sumstats`),
#'   'uniform' (uniform prior),
#'   or 'custom' (use the values provided in \code{prior_weights}).
#' @param prior_weights A vector of prior probability for each SNP.
#' @param L Number of causal signals.
#' If L = 1, bigSNP or region_info are not required.
#' @param estimate_residual_variance The default is FALSE,
#' the residual variance is fixed to 1 or variance of y.
#' If the in-sample LD matrix is provided,
#' we recommend setting \code{estimate_residual_variance = TRUE}.
#' @param verbose If TRUE, print progress.
#' @param save If TRUE, save SuSiE result and LD (R) matrix for each locus.
#' @param outputdir Directory of SuSiE result
#' @param outname Filename of SuSiE result
#' @return A list of SuSiE results; one per LD block.
#' @export
run_finemapping <- function(sumstats,
                            bigSNP = NULL,
                            region_info = NULL,
                            n = NULL,
                            priortype = c('torus', 'uniform', 'custom'),
                            prior_weights = NULL,
                            L = 1,
                            estimate_residual_variance = FALSE,
                            verbose = FALSE,
                            save = FALSE,
                            outputdir = getwd(),
                            outname = NULL,
                            ...){

  priortype <- match.arg(priortype)

  if(priortype == 'torus'){
    stopifnot('torus_prior' %in% colnames(sumstats))
    prior_weights <- sumstats$torus_prior
  }else if(priortype == 'uniform'){
    prior_weights <- NULL
  }else if(priortype == 'custom'){
    cat('use custom prior_weights \n')
    stopifnot(length(prior_weights) == nrow(sumstats))
  }

  finemap.locus.list <- unique(sumstats$locus)
  cat(sprintf("Finemapping %d loci...\n", length(finemap.locus.list)))

  susie_results <- list()
  for(locus in finemap.locus.list){
    cat(sprintf('Finemapping locus %s...\n', locus))
    sumstats_locus <- sumstats[sumstats$locus == locus, ]

    if (!is.null(region_info)) {
      # load precomputed LD (R) matrix and variant info
      LD_matrix_file <- region_info$LD_matrix[region_info$locus == locus]
      if (!file.exists(LD_matrix_file))
        stop(paste('LD matrix file', LD_matrix_file, 'does not exist!'))

      snp_info_file <- region_info$snp_info[region_info$locus == locus]
      if (!file.exists(snp_info_file))
        stop(paste('LD SNP info file', snp_info_file, 'does not exist!'))

      R <- read_LD(LD_matrix_file)
      snp_info <- data.table::fread(snp_info_file)

      res <- susie_finemap_region(sumstats_locus,
                                  R=R,
                                  snp_info = snp_info,
                                  n=n,
                                  L=L,
                                  prior_weights=prior_weights,
                                  estimate_residual_variance=estimate_residual_variance,
                                  verbose=verbose,
                                  save = save,
                                  ...)

    } else if (!is.null(bigSNP)) {
      # compute LD from bigSNP genotype panel
      res  <- susie_finemap_region(sumstats_locus,
                                   bigSNP=bigSNP,
                                   n=n,
                                   L=L,
                                   prior_weights=prior_weights,
                                   estimate_residual_variance=estimate_residual_variance,
                                   verbose=verbose,
                                   save = save,
                                   ...)

    } else if (L == 1) {
      res <- susie_finemap_region(sumstats_locus,
                                  n=n,
                                  L=1,
                                  prior_weights=prior_weights,
                                  estimate_residual_variance=estimate_residual_variance,
                                  verbose=verbose,
                                  save = FALSE,
                                  ...)
    } else {
      stop("Please provide region_info or bigSNP object when L > 1")
    }

    susie_results[[as.character(locus)]] <- res

    if (isTRUE(save)) {
      if (!dir.exists(outputdir))
        dir.create(outputdir, recursive = TRUE)

      saveRDS(res, file = file.path(outputdir, paste0(outname, '.locus', locus, '.susie.res.rds')))
      rm(res)
    }

    if(verbose){
      cat(sprintf('%.0f%% completed.\n', length(susie_results)/length(finemap.locus.list)*100))
    }

  }

  return(susie_results)

}

#' @title merges SuSiE results with original summary statistics data frame
#' @description  merges SuSiE results with original summary statistics data frame
#' @param susie_results SuSiE finemapping result
#' @param sumstats A data frame of summary statistics
#' @return A data frame with summary statistics and two additional columns
#' of 'susie_pip' and 'cs' (CS indices)
#' @export
merge_susie_sumstats <- function(susie_results, sumstats){

  sumstats$susie_pip <- NA
  sumstats$cs <- NA

  if (class(susie_results) == "susie"){
    susie_results <- list(susie_results)
  }

  for(i in 1:length(susie_results)){

    susie.locus.res <- susie_results[[i]]
    stopifnot(class(susie.locus.res) == "susie")

    # add PIP result of the SNPs in this locus
    n.snps <- length(susie.locus.res$pip)
    snp.idx <- match(susie.locus.res$id, sumstats$snp)
    sumstats[snp.idx, 'susie_pip'] <- susie.locus.res$pip

    # add CS index of the SNPs in this locus
    # set cs.index = 0 for SNPs not in credible sets
    cs.index <- rep(0, n.snps)
    snp.idx.inCS <- unlist(susie.locus.res$sets$cs)
    if( length(snp.idx.inCS) > 0 ){
      cs.index[snp.idx.inCS] <- rep(susie.locus.res$sets$cs_index,
                                    lapply(susie.locus.res$sets$cs, length))
    }
    sumstats[snp.idx, 'cs'] <-  cs.index
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

