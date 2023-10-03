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
                                  LD_Blocks,
                                  bigSNP){

  cat('Cleaning summary statistics...\n')
  cleaned.sumstats <- clean_sumstats(sumstats,
                                     chr=chr, pos=pos, beta=beta, se=se,
                                     a0=a0, a1=a1, snp=snp, pval=pval)

  if(missing(LD_Blocks)){
    cat('Skipped assigning SNPs to LD blocks. \n')
  }else{
    cat('Assigning GWAS SNPs to LD blocks...\n')
    cleaned.sumstats <- assign_snp_locus(cleaned.sumstats, LD_Blocks)
  }

  if(missing(bigSNP)){
    cat('Skipped matching GWAS with bigSNP reference panel. \n')
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

  # Check chromosome names, remove 'chr'
  if( any(grepl('chr', cleaned.sumstats$chr)) ){
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


#' @title Assign GWAS SNPs to LD blocks
#' @param sumstats A data frame of GWAS summary statistics.
#' @param LD_Blocks A data frame of LD blocks with four columns,
#' 'chr', 'start', 'end', and 'locus'.
#' @export
assign_snp_locus <- function(sumstats, LD_Blocks){

  LD_Blocks <- as.data.frame(LD_Blocks)
  colnames(LD_Blocks)[1:4] <- c('chr', 'start', 'end', 'locus')

  if( any(grepl('chr', sumstats$chr)) ){
    sumstats$chr <- gsub('chr', '', sumstats$chr)
  }

  if( any(grepl('chr', LD_Blocks$chr)) ){
    LD_Blocks$chr <- gsub('chr', '', LD_Blocks$chr)
  }

  LD_Blocks.gr <- GenomicRanges::makeGRangesFromDataFrame(LD_Blocks,
                                                          keep.extra.columns = TRUE)

  snp.gr <- GenomicRanges::makeGRangesFromDataFrame(sumstats,
                                                    start.field = 'pos', end.field = 'pos')

  snp.gr <- plyranges::mutate(snp.gr, snp=sumstats$snp)

  snp.ld.block.overlap <- plyranges::join_overlap_inner(snp.gr, LD_Blocks.gr)
  snp.ld.block <- tibble::as_tibble(snp.ld.block.overlap@elementMetadata)
  # remove duplicated SNPs,
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ]
  sumstats.ld.block <- dplyr::inner_join(sumstats, snp.ld.block, by = 'snp')

  return(sumstats.ld.block)
}

#' @title Match GWAS with bigSNP reference panel.
#' @param sumstats  A data frame of GWAS summary statistics
#' @param bigSNP  bigSNP object
#' @param strand_flip Whether to try to flip strand? (default is TRUE).
#' If so, ambiguous alleles A/T and C/G are removed.
#' @param match.min.prop Minimum proportion of variants in the smallest data
#' to be matched, otherwise stops with an error. Default: 10%
#' @export
match_gwas_bigsnp <- function(sumstats, bigSNP, strand_flip = TRUE, match.min.prop = 0.1){

  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')

  matched.sumstats <- bigsnpr::snp_match(sumstats,
                                         snp_info,
                                         strand_flip = strand_flip,
                                         match.min.prop = match.min.prop)

  matched.sumstats <- matched.sumstats %>%
    tibble::as_tibble() %>%
    dplyr::rename(og_index = `_NUM_ID_.ss`) %>%
    dplyr::rename(bigSNP_index = `_NUM_ID_`) %>%
    dplyr::mutate(zscore = beta/se)

  return(matched.sumstats)
}
