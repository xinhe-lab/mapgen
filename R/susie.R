
#' @title PrepareSusieData
#' @description Adds torus results to cleaned summary statistics
#' @param sumstats a tibble or data frame containing raw summary statistics
#' @param torus_pip a tibble containing PIP of each SNP (result from RunTorus)
#' @param torus_fdr a tibble containing the FDR of each region (result from RunTorusFDR)
#' @return tibble of summary statistics updated with torus output
#' @export
PrepareSusieData <- function(sumstats, torus_pip, torus_fdr, fdr_thresh=0.1){

  # keep loci at fdr_thresh FDR (10% by default)
  chunks <- torus_fdr$region_id[torus_fdr$fdr < fdr_thresh]
  sumstats <- sumstats[sumstats$locus %in% chunks, ]

  # Add Torus PIP
  sumstats <- dplyr::inner_join(sumstats, torus_pip, by='snp')

  return(sumstats)

}

#' @title RunFinemapping
#' @description Runs SuSiE with L = 1
#' @param sumstats a tibble or data frame containing raw summary statistics; must have header!
#' @param bigSNP a bigsnpr object attached via bigsnpr::snp_attach()
#' @param priortype prior type: "torus" or "uniform".
#' @return list of finemapping results; one per LD block
#' @export
RunFinemapping <- function(sumstats, bigSNP, priortype = c('torus', 'uniform')){

  stopifnot('torus_pip' %in% colnames(sumstats))
  priortype <- match.arg(priortype)

  if(priortype == 'torus'){
    usePrior <- TRUE
  }else if(priortype == 'uniform'){
    usePrior <- FALSE
  }

  chunks <- unique(sumstats$locus)

  susie_res <- list()
  for(z in seq_along(chunks)){
    cat(sprintf('Finemapping chunks %d of %d ...\n', z, length(chunks)))
    susie.df <- sumstats[sumstats$locus == z, ]
    susie_res[[as.character(z)]] <- run.susie(susie.df, bigSNP, z, L = 1, prior = usePrior)
  }

  return(susie_res)

}


#' @title run SUSIE
#' @param sumstats summary statistics
#'
#' @param bigSNP bigSNP object
#' @param ldchunk LD chunk
#' @param L Number of causal signals
#' @param prior Logical, if TRUE, use the \code{torus_pip} column
#' in \code{sumstats} as prior
#'
#' @export
run.susie <- function(sumstats, bigSNP, ldchunk, L, prior){

  sub.sumstats <- sumstats[sumstats$locus == ldchunk, ]
  if(nrow(sub.sumstats) > 1){
    X <- bigSNP$genotypes[ , sub.sumstats$bigSNP_index]
    X <- scale(X, center = T, scale = T)
    zhat <- sub.sumstats$zscore
    R <- cov2cor((crossprod(X) + tcrossprod(zhat))/nrow(X))
    if(prior){
      res <- suppressWarnings(susieR::susie_rss(z = zhat,
                                                prior_weights = sub.sumstats$torus_pip,
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


