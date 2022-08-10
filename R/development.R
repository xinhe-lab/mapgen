
#' @title Run fine-mapping with SuSiE for one LD block
#' @param sumstats summary statistics
#' @param R p x p LD (correlation) matrix
#' @param bigSNP bigSNP object
#' @param locus LD block index
#' @param L Number of causal signals
#' @param useprior Logical, if TRUE, use the \code{torus_prior} column
#' in \code{sumstats} as prior.
#' @param estimate_residual_variance The default is FALSE,
#' the residual variance is fixed to 1 or variance of y.
#' If the in-sample LD matrix is provided,
#' we recommend setting estimate_residual_variance = TRUE.
#' @return finemapping results for one LD block
run_susie_v2 <- function(sumstats, R, bigSNP, locus, L=1, useprior,
                      estimate_residual_variance=FALSE){

  sub.sumstats <- sumstats[sumstats$locus == locus, ]
  if(nrow(sub.sumstats) > 1){
    if(missing(R) || is.null(R)){
      if(missing(bigSNP) || is.null(bigSNP)){
        stop("Please provide R matrix or bigSNP object!")
      }
      cat("Computing R from bigSNP genotype matrix...\n")
      X <- bigSNP$genotypes[ , sub.sumstats$bigSNP_index]
      X <- scale(X, center = T, scale = T)
      R <- cor(X)
    }
    z <- sub.sumstats$zscore

    if(useprior){
      prior_weights <- sub.sumstats$torus_prior
    }else{
      prior_weights <- NULL
    }

    res <- suppressWarnings(susieR::susie_rss(z = z,
                                              prior_weights = prior_weights,
                                              R = R,
                                              L = L,
                                              estimate_residual_variance = estimate_residual_variance,
                                              verbose = F))
    return(res)
  }
}

# Run fine-mapping with SuSiE with reference LD matrix
run_susie_LD <- function(sumstats, R, locus, L=1, useprior,
                         estimate_residual_variance=FALSE){

  sub.sumstats <- sumstats[sumstats$locus == locus, ]
  if(nrow(sub.sumstats) > 1){
    z <- sub.sumstats$zscore
    if(useprior){
      prior_weights <- sub.sumstats$torus_prior
    }else{
      prior_weights <- NULL
    }
    res <- suppressWarnings(susieR::susie_rss(z = z,
                                              prior_weights = prior_weights,
                                              R = R,
                                              L = L,
                                              estimate_residual_variance = estimate_residual_variance,
                                              verbose = F))
    return(res)
  }
}
