# check allele flipping
# Copied from allele.qc function
# in https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R,
allele.qc = function(a1,a2,ref1,ref2) {

  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

  return(snp)
}

#' Harmonize GWAS summary statistics with LD reference
#'
#' @param sumstats A data frame of GWAS summary statistics,
#' including columns "snp" (SNP ID), a0" (reference allele), "a1" (effect allele), "beta", and "zscore".
#'
#' @param LD_snp_info a data frame, SNP info for LD reference,
#'  with columns "chrom", "id", "pos", "alt", "ref", "locus" (optional)
#'
#' @param strand_flip Whether to flip signs when reverse complement matches? (default is TRUE).
#'
#' @param remove_strand_ambig Whether to remove ambiguous alleles (A/T and C/G)? (default is TRUE).
#'
#' @return a data frame of harmonized GWAS summary statistics
#'
harmonize_sumstats_LD <- function(sumstats,
                                  LD_snp_info,
                                  strand_flip = TRUE,
                                  remove_strand_ambig = TRUE){

  cat(sprintf("%d variants in sumstats before harmonization...\n", nrow(sumstats)))

  # keep GWAS SNPs in LD reference
  sumstats <- sumstats[sumstats$snp %in% LD_snp_info$id,]

  if (length(sumstats$snp) == 0)
    stop("No SNPs found in both sumstats and LD reference")

  cat(sprintf("%d variants in both sumstats and LD reference ...\n", length(sumstats$snp)))

  LD.idx <- match(sumstats$snp, LD_snp_info$id)
  LD_snp_info <- LD_snp_info[LD.idx, ]

  if (!is.null(LD_snp_info$locus)){
    sumstats$locus <- LD_snp_info$locus
  }

  # QC / allele-flip the input and output
  qc <- allele.qc(sumstats$a1, sumstats$a0, LD_snp_info$alt, LD_snp_info$ref)
  flip.idx <- which(qc[["flip"]])
  remove.idx <- which(!qc[["keep"]])

  # Flip signs for allele flipped variants
  if (isTRUE(strand_flip) & any(flip.idx)) {
    cat(sprintf("Flip signs for %d (%.2f%%) variants with allele flipped. \n",
                length(flip.idx), length(flip.idx)/length(sumstats$snp)*100))

    sumstats[flip.idx, c("a0", "a1")] <- sumstats[flip.idx, c("a1", "a0")]
    if("beta" %in% colnames(sumstats)){
      sumstats[flip.idx, "beta"] <- -sumstats[flip.idx, "beta"]
    }
    if("zscore" %in% colnames(sumstats)){
      sumstats[flip.idx, "zscore"] <- -sumstats[flip.idx, "zscore"]
    }else{
      sumstats$zscore <- sumstats$beta / sumstats$se
    }
  }

  # Remove strand ambiguous SNPs (if any)
  if (isTRUE(remove_strand_ambig) && any(remove.idx)) {
    sumstats <- sumstats[-remove.idx, , drop = F]
    cat(sprintf("Remove %d (%.2f%%) strand ambiguous variants. \n",
                length(remove.idx), length(remove.idx)/length(sumstats$snp)*100))
  }

  return(sumstats)
}


#' @title Match alleles between GWAS summary statistics and bigSNP reference panel.
#' @description
#' Match alleles between summary statistics and SNP information in the
#' bigSNP reference panel using the \code{bigsnpr::snp_match()} function.
#' Match by ("chr", "a0", "a1") and ("pos" or "rsid"),
#' accounting for possible strand flips and reverse reference alleles (opposite effects).
#'
#' @param sumstats  A data frame of GWAS summary statistics,
#' with columns "chr", "pos", "a0", "a1" and "beta".
#' @param bigSNP a \code{bigsnpr} object attached via \code{bigsnpr::snp_attach()}
#' containing the reference genotype panel.
#' @param strand_flip Whether to try to flip strand? (default is TRUE).
#' If so, ambiguous alleles A/T and C/G are removed.
#' @param match.min.prop Minimum proportion of variants in the smallest data
#' to be matched, otherwise stops with an error. Default: 10%
#' @return A data frame with matched summary statistics.
#' Values in column "beta" are multiplied by -1 for variants with
#' alleles reversed (i.e. swapped).
#' New variable "ss_index" returns the corresponding row indices of the sumstats,
#' and "bigSNP_index" corresponding to the indices of the bigSNP.
#' @export
match_gwas_bigsnp <- function(sumstats,
                              bigSNP,
                              strand_flip = TRUE,
                              match.min.prop = 0.1, ...){

  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')

  matched.sumstats <- bigsnpr::snp_match(sumstats,
                                         snp_info,
                                         strand_flip = strand_flip,
                                         match.min.prop = match.min.prop,
                                         ...)

  matched.sumstats <- matched.sumstats %>%
    tibble::as_tibble() %>%
    dplyr::rename(ss_index = `_NUM_ID_.ss`) %>%
    dplyr::rename(bigSNP_index = `_NUM_ID_`) %>%
    dplyr::mutate(zscore = beta/se)

  return(matched.sumstats)
}
