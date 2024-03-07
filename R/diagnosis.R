
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
#' @export
LD_diagnosis_susie_rss <- function(z, R, n){
  cat('Estimate consistency between the z-scores and reference LD matrix ...\n')
  lambda <- susieR::estimate_s_rss(z = z, R = R, n = n)
  cat('Estimated lambda =', lambda, '\n')

  cat('Compute expected z-scores based on conditional distribution of z-scores ...\n')
  condz <- susieR::kriging_rss(z = z, R = R, n = n, s = lambda)

  return(condz)
}
