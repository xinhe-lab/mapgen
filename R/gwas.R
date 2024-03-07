#' @title Cleans GWAS summary statistics and adds metadata
#' @description Cleans GWAS summary statistics and adds metadata, including:
#' the index in the bigSNP object of each SNP, and the LD block locus.
#'
#' @param sumstats A data frame of GWAS summary statistics.
#' @param chr Name of the chromosome column in summary statistics.
#' @param pos Name of the position column (base pair position).
#' @param beta Name of beta column (if you have Odds Ratio,
#' you will need to transform it to log(Odds Ratio)).
#' @param se Name of the standard error (se) column.
#' @param a0 Column name of the reference allele.
#' @param a1 Column name of the association/effect allele.
#' @param snp Name of the SNP ID (rsID) column.
#' @param pval Name of the p-value column.
#' @param remove_indels If TRUE, remove indels
#' @param bigSNP a \code{bigsnpr} object attached via \code{bigsnpr::snp_attach()}
#' containing the reference genotype panel.
#' @param LD_Blocks A data frame of LD blocks with four columns,
#' 'chr', 'start', 'end', and 'locus' (LD block indices).
#' @return A data frame of processed GWAS summary statistics.
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
                                  remove_indels = TRUE,
                                  LD_Blocks = NULL,
                                  bigSNP = NULL,
                                  region_info = NULL,
                                  strand_flip = TRUE,
                                  remove_strand_ambig = TRUE,
                                  ...){

  cat('Cleaning summary statistics...\n')
  sumstats <- clean_sumstats(sumstats,
                             chr=chr, pos=pos, beta=beta, se=se,
                             a0=a0, a1=a1, snp=snp, pval=pval,
                             remove_indels=remove_indels)

  if(!missing(LD_Blocks)){
    cat('Assigning GWAS SNPs to LD blocks...\n')
    sumstats <- assign_snp_locus(sumstats, LD_Blocks)
  }

  if(!missing(bigSNP)){
    cat('Matching GWAS with bigSNP reference panel...\n')
    sumstats <- match_gwas_bigsnp(sumstats,
                                  bigSNP,
                                  strand_flip = strand_flip,
                                  ...)

  } else if(!missing(region_info)){
    cat('Harmonizing GWAS with LD reference panel...\n')

    # Read SNP info for all LD regions
    cat('Reading LD reference SNPs ... \n')
    LD_snp_info.list <- lapply(1:nrow(region_info), function(i){
      df <- data.table::fread(region_info$snp_info[i])
      df$locus <- region_info$locus[i]
      df
    })
    LD_snp_info <- do.call(rbind, LD_snp_info.list)

    sumstats <- harmonize_sumstats_LD(sumstats,
                                      LD_snp_info,
                                      strand_flip = strand_flip,
                                      remove_strand_ambig = remove_strand_ambig)

  }

  return(sumstats)

}

#' @title Cleans summary statistics
#' @description
#' It will extract the required columns from summary statistics,
#' check chromosomes, remove X, Y chromosomes, compute z-scores,
#' convert alleles to upper case, remove indels,
#' and sort by chromosome and position.
#' @param sumstats A data frame of GWAS summary statistics.
#' It is required to have the following columns:
#' chr, position, beta, se, a0, a1, SNP ID (rs), p-value.
#' @param chr Name of the chromosome column in summary statistics.
#' @param pos Name of the position column (base pair position).
#' @param beta Name of beta column (if you have Odds Ratio,
#' you will need to transform it to log(Odds Ratio)).
#' @param se Name of the standard error (se) column.
#' @param a0 Column name of the reference allele.
#' @param a1 Column name of the association/effect allele.
#' @param snp Name of the SNP ID (rsID) column.
#' @param pval Name of the p-value column.
#' @return A data frame of cleaned summary statistics,
#' sort by chromosome and position.
#' @export
clean_sumstats <- function(sumstats,
                           chr = 'chr',
                           pos = 'pos',
                           beta = 'beta',
                           se = 'se',
                           a0 = 'a0',
                           a1 = 'a1',
                           snp = 'snp',
                           pval = 'pval',
                           remove_indels = TRUE){

  cols.to.keep <- c(chr, pos, beta, se, a0, a1, snp, pval)

  if(!all(cols.to.keep %in% colnames(sumstats))){
    stop("sumstats needs to contain the following columns: ",
         paste(cols.to.keep, collapse = " "))
  }

  # Extract relevant columns
  cleaned.sumstats <- sumstats[, cols.to.keep]
  colnames(cleaned.sumstats) <- c('chr','pos','beta','se','a0','a1','snp','pval')

  # Check chromosomes
  # Remove 'chr'
  if( any(grepl('chr', cleaned.sumstats$chr)) ){
    cleaned.sumstats$chr <- gsub('chr', '', cleaned.sumstats$chr)
  }

  # Remove X, Y chromosomes
  cleaned.sumstats <- cleaned.sumstats[!(cleaned.sumstats$chr %in% c('X','Y')), ]
  cleaned.sumstats$chr <- as.integer(cleaned.sumstats$chr)

  # Compute z-scores
  zscore <- cleaned.sumstats$beta/cleaned.sumstats$se
  cleaned.sumstats$zscore <- zscore
  cleaned.sumstats <- cleaned.sumstats[!is.na(zscore),]

  # Convert alleles to upper case
  cleaned.sumstats$a0 <- toupper(cleaned.sumstats$a0)
  cleaned.sumstats$a1 <- toupper(cleaned.sumstats$a1)

  # Keep SNPs only, remove indels
  if(remove_indels){
    nucs <- c('A','C','T','G')
    cleaned.sumstats <- cleaned.sumstats %>% dplyr::filter(a0 %in% nucs, a1 %in% nucs)
  }

  # Sort by chromosome and position
  cleaned.sumstats <- cleaned.sumstats %>% dplyr::arrange(chr, pos)

  # Remove duplicate SNPs
  chr_pos <- paste0(cleaned.sumstats$chr, '_', cleaned.sumstats$pos)
  cleaned.sumstats <- cleaned.sumstats[!duplicated(chr_pos), ]

  return(cleaned.sumstats)
}


#' @title Assign GWAS SNPs to LD blocks
#' @param sumstats A data frame of GWAS summary statistics.
#' It is required to have the following columns:
#' chr, pos, snp (rsID).
#' @param LD_Blocks A data frame of LD blocks with four columns,
#' 'chr', 'start', 'end', and 'locus' (LD block indices).
#' @return A data frame with summary statistics with assigned locus ID.
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

#' Match GWAS sumstats with LD reference files. Only keep variants included in
#' LD reference.
#'
#' @param sumstats A data frame of GWAS summary statistics.
#' @param R LD matrix
#' @param snp_info Variant information for the LD matrix.
#'
#' @return A list, containing matched GWAS summary statistics and LD matrix.
#' @export
match_gwas_LDREF <- function(sumstats, R, snp_info){
  stopifnot(nrow(R) == nrow(snp_info))
  sumstats <- sumstats[sumstats$snp %in% snp_info$id,]
  LD.idx <- match(sumstats$snp, snp_info$id)
  R <- R[LD.idx, LD.idx]
  snp_info <- snp_info[LD.idx, ]
  return(list(sumstats = sumstats, R = R, snp_info = snp_info))
}
