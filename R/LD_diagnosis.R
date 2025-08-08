
#' Perform diagnosis to check the consistency between the GWAS z-scores and LD matrix,
#' and detect problematic z-scores (e.g. allele switch issue) using
#' diagnostic functions from SuSiE RSS.
#'
#' @param z z scores.
#' @param R LD correlation matrix.
#' @param n The sample size. (Optional, but highly recommended.)
#'
#' @return a data frame of expected z-scores and test statistics
#' from \code{susieR::kriging_rss}.
#' @importFrom susieR estimate_s_rss kriging_rss
#' @export
LD_diagnosis_susie_rss <- function(z, R, n){

  cat('Estimate consistency between the z-scores and reference LD matrix ...\n')
  lambda <- susieR::estimate_s_rss(z = z, R = R, n = n)
  cat('Estimated lambda =', lambda, '\n')

  # Compute expected z-scores based on conditional distribution of z-scores
  cat('Compute expected z-scores based on conditional distribution of z-scores ...\n')
  kriging_res <- susieR::kriging_rss(z = z, R = R, n = n, s = lambda)

  # compute p-values for the significance of z-score difference between observed and estimated values
  kriging_res$conditional_dist$p_diff <- pchisq(kriging_res$conditional_dist$z_std_diff^2, df = 1, lower.tail=FALSE)

  return(kriging_res)
}
