#' @title Cleans GWAS summary statistics and adds metadata
#' @description Cleans GWAS summary statistics and adds metadata, including:
#' the locus ID of every SNP, as well as its index in the bigSNP object.
#'
#' @param sumstats A data frame of GWAS summary statistics.
#' @param chr Name of the chromosome column in summary statistics file.
#' @param pos Name of the position column (base pair position).
#' @param beta Name of beta column (if you have Odds Ratio,
#' you will need to transform it to log(Odds Ratio)).
#' @param se Name of the standard error (se) column.
#' @param a0 Column name of the reference allele.
#' @param a1 Column name of the association/effect allele.
#' @param snp Name of the SNP ID (rsID) column.
#' @param pval Name of the p-value column.
#' @param bigSNP bigSNP object from \code{bigsnpr} containing the reference
#' genotype panel.
#' @param LD_Blocks Reference LD blocks.
#'
#' @return Processed GWAS summary statistics.
#' @export
process_gwas_sumstats <- function(sumstats,
                                  chr = 'chr',
                                  pos = 'pos',
                                  beta = 'beta',
                                  se = 'se',
                                  a0 = 'a0',
                                  a1 = 'a1',
                                  snp = 'snp',
                                  pval = 'pval',
                                  LD_Blocks = NULL,
                                  bigSNP = NULL){

  cat('Cleaning summary statistics...\n')
  cleaned.sumstats <- clean_sumstats(sumstats,
                                     chr=chr, pos=pos, beta=beta, se=se,
                                     a0=a0, a1=a1, snp=snp, pval=pval)

  if(is.null(LD_Blocks)){
    cat('LD_Blocks not provided. Skipped assigning SNPs to LD blocks. \n')
  }else{
    cat('Assigning GWAS SNPs to LD blocks...\n')
    cleaned.sumstats <- assign_snp_locus(cleaned.sumstats, LD_Blocks)
  }

  if(is.null(bigSNP)){
    cat('bigSNP not provided. Skipped matching GWAS with bigSNP reference panel. \n')
  }else{
    cat('Matching GWAS with bigSNP reference panel...\n')
    cleaned.sumstats <- match_gwas_bigsnp(cleaned.sumstats, bigSNP)
  }

  return(cleaned.sumstats)

}

#' @title Cleans summary statistics
#' @param sumstats A data frame of GWAS summary statistics.
#' @param cols.to.keep columns to keep.
#' It is required to have 8 columns in the exact order of:
#' chr, position, beta, se, allele1, allele2, SNP ID (rs), p-value.
#'
#' @export
clean_sumstats <- function(sumstats,
                           chr = 'chr',
                           pos = 'pos',
                           beta = 'beta',
                           se = 'se',
                           a0 = 'a0',
                           a1 = 'a1',
                           snp = 'snp',
                           pval = 'pval'){

  stopifnot(!is.null(sumstats))

  cols.to.keep <- c(chr, pos, beta, se, a0, a1, snp, pval)

  if(!all(cols.to.keep %in% colnames(sumstats))){
    stop(sprintf('Column: %s cannot be found in the summary statistics!',
                 cols.to.keep[which(!cols.to.keep %in% colnames(sumstats))]))
  }else{
    # Extract relevant columns
    cleaned.sumstats <- sumstats[, cols.to.keep]
    colnames(cleaned.sumstats) <- c('chr','pos','beta','se','a0','a1','snp','pval')
  }

  # Check chromosome names
  if( any(grepl('chr', cleaned.sumstats$chr)) ){
    cat("Remove \'chr\' from the chr column...\n")
    cleaned.sumstats$chr <- gsub('chr', '', cleaned.sumstats$chr)
  }

  # Remove X, Y chromosomes
  cleaned.sumstats <- cleaned.sumstats[!(cleaned.sumstats$chr %in% c('X','Y')), ]
  cleaned.sumstats$chr <- as.integer(cleaned.sumstats$chr)

  # Compute Zscores
  zscore <- cleaned.sumstats$beta/cleaned.sumstats$se
  cleaned.sumstats$zscore <- zscore
  cleaned.sumstats <- cleaned.sumstats[!is.na(zscore),]

  # Convert alleles to upper case
  cleaned.sumstats$a0 <- toupper(cleaned.sumstats$a0)
  cleaned.sumstats$a1 <- toupper(cleaned.sumstats$a1)

  # Keep SNPs only, no indels
  nucs <- c('A','C','T','G')
  bola1 <- (cleaned.sumstats$a0 %in% nucs)
  bola2 <- (cleaned.sumstats$a1 %in% nucs)
  cleaned.sumstats <- cleaned.sumstats[bola1 & bola2,]

  # Sort by chromosome and position
  cleaned.sumstats <- cleaned.sumstats[order(cleaned.sumstats$chr, cleaned.sumstats$pos), ]

  # Remove duplicate SNPs
  chrpos <- paste0(cleaned.sumstats$chr, '_', cleaned.sumstats$pos)
  cleaned.sumstats <- cleaned.sumstats[!duplicated(chrpos), ]

  return(cleaned.sumstats)
}


#' @title Assigns each SNP to one LD block
#' @param cleaned.sumstats A data frame of GWAS summary statistics.
#' @param LD_Blocks A data frame of LD blocks
#' @export
assign_snp_locus <- function(cleaned.sumstats, LD_Blocks){

  ld.gr <- make_ranges(LD_Blocks$X1, LD_Blocks$X2, LD_Blocks$X3)
  ld.gr <- plyranges::mutate(ld.gr, locus=LD_Blocks$X4)

  snp.gr <- make_ranges(seqname = cleaned.sumstats$chr,
                           start = cleaned.sumstats$pos,
                           end = cleaned.sumstats$pos)
  snp.gr <- plyranges::mutate(snp.gr, snp=cleaned.sumstats$snp)

  snp.ld.overlap <- plyranges::join_overlap_inner(snp.gr, ld.gr)
  snp.ld.block <- tibble::as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple LD blocks due to edge of LD blocks
  cleaned.annot.sumstats <- dplyr::inner_join(cleaned.sumstats, snp.ld.block, 'snp')

  return(cleaned.annot.sumstats)
}

#' @title Match GWAS with bigSNP reference panel.
#' @param sumstats  A data frame of GWAS summary statistics
#' @param bigSNP  bigSNP object
#' @param strand_flip Whether to try to flip strand? (default is TRUE).
#' If so, ambiguous alleles A/T and C/G are removed.
#' @param match.min.prop Minimum proportion of variants in the smallest data
#' to be matched, otherwise stops with an error. Default: 10%
#' @importFrom magrittr %>%
#' @export
match_gwas_bigsnp <- function(sumstats, bigSNP, strand_flip = TRUE, match.min.prop = 0.1){

  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')

  matched.sumstats <- tibble::as_tibble(bigsnpr::snp_match(sumstats,
                                                           snp_info,
                                                           strand_flip = strand_flip,
                                                           match.min.prop = match.min.prop))

  matched.sumstats <- matched.sumstats %>%
    dplyr::rename(og_index = `_NUM_ID_.ss`) %>%
    dplyr::rename(bigSNP_index = `_NUM_ID_`) %>%
    dplyr::mutate(zscore = beta/se)

  return(matched.sumstats)
}
