#' @title Cleans GWAS summary statistics and adds metadata
#' @description Cleans GWAS summary statistics and adds metadata
#'
#' @param sumstats_file file containing raw summary statistics.
#' Coordinates should be hg19/b37; we current do not support hg38.
#' @param chr Name of the chromosome column in summary statistics file.
#' @param pos Name of the position column (base pair position, in *hg19*!).
#' @param beta Name of beta column (if you have Odds Ratio,
#' you will need to transform it to log(Odds Ratio)).
#' @param se Name of the standard error (se) column.
#' @param a0 Column name of the reference allele.
#' @param a1 Column name of the association/effect allele.
#' @param snp Name of the SNP ID (rsID) column.
#' @param pval Name of the p-value column.
#' @param bigSNP bigSNP object from \code{bigsnpr}.
#' @param LD_Blocks Reference LD blocks
#'
#' @return Cleaned summary statistics + LD block of every SNP,
#' as well as its index in the reference panel of genotypes
#' @export
process_gwas_sumstats <- function(sumstats_file,
                                  chr = 'chr',
                                  pos = 'pos',
                                  beta = 'beta',
                                  se = 'se',
                                  a0 = 'a0',
                                  a1 = 'a1',
                                  snp = 'snp',
                                  pval = 'pval',
                                  bigSNP,
                                  LD_Blocks){

  cat('Loading summary statistics...\n')
  if(is.character(sumstats_file) & length(sumstats_file) == 1){
    sumstats <- vroom::vroom(sumstats_file, col_names = TRUE, show_col_types = FALSE)
  }

  cat('Cleaning summary statistics...\n')
  cleaned_sumstats <- clean_sumstats(sumstats,
                                     chr=chr, pos=pos, beta=beta, se=se,
                                     a0=a0, a1=a1, snp=snp, pval=pval)

  cat('Matching to reference panel...\n')
  cleaned_sumstats <- merge_bigsnp_gwas(cleaned_sumstats, bigSNP = bigSNP)

  cat('Assigning SNPs to LD blocks...\n')
  if(missing(LD_Blocks)){
    cat('No LD blocks supplied. Using the included European LD blocks from 1KG....\n')
    data('Euro_LD_Chunks', package='mapgen')
  }
  cleaned_sumstats <- assign_locus_snp(cleaned.sumstats = cleaned_sumstats,
                                       ld = LD_Blocks)

  return(cleaned_sumstats)

}

#' @title Cleans summary statistics
#' @param sumstats a tibble or data frame containing raw summary statistics.
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
    cat('Rename and extract columns in summary statistics ...\n')
    cleaned.sumstats <- sumstats[, cols.to.keep]
    colnames(cleaned.sumstats) <- c('chr','pos','beta','se','a0','a1','snp','pval')
  }

  # Check chromosome names
  if( any(grepl('chr', cleaned.sumstats$chr)) ){
    cat("Remove \'chr\' from the chr column...\n")
    cleaned.sumstats$chr <- gsub('chr', '', cleaned.sumstats$chr)
  }

  # drop X, Y chromosomes
  cat('Drop X,Y chromosomes ...\n')
  cleaned.sumstats <- cleaned.sumstats[!(cleaned.sumstats$chr %in% c('X','Y')), ]
  # make chromosomes integers
  cat('Chromosomes included: ', unique(cleaned.sumstats$chr), '\n')
  cleaned.sumstats$chr <- as.integer(cleaned.sumstats$chr)

  # Compute Zscores
  cat('Compute z-scores ...\n')
  zscore <- cleaned.sumstats$beta/cleaned.sumstats$se
  cleaned.sumstats['zscore'] <- zscore
  cleaned.sumstats <- cleaned.sumstats[!is.na(zscore),]

  # convert alleles to upper case
  cleaned.sumstats$a0 <- toupper(cleaned.sumstats$a0)
  cleaned.sumstats$a1 <- toupper(cleaned.sumstats$a1)

  # Keep SNPs only, no indels
  cat('Remove indels ...\n')
  nucs <- c('A','C','T','G')
  bola1 <- (cleaned.sumstats$a0 %in% nucs)
  bola2 <- (cleaned.sumstats$a1 %in% nucs)
  cleaned.sumstats <- cleaned.sumstats[bola1 & bola2,]

  # Sort by chromosome and position
  cleaned.sumstats <- cleaned.sumstats[order(cleaned.sumstats$chr, cleaned.sumstats$pos), ]

  # drop duplicate SNPs
  cat('Remove duplicate SNPs ...\n')
  chrpos <- paste0(cleaned.sumstats$chr, '_', cleaned.sumstats$pos)
  cleaned.sumstats <- cleaned.sumstats[!duplicated(chrpos), ]

  return(cleaned.sumstats)
}


#' @title Assigns each SNP to one ld-block
#' @param cleaned.sumstats Cleaned summary statistics
#' @param ld LD blocks
#' @importFrom magrittr %>%
#' @export
assign_locus_snp <- function(cleaned.sumstats, ld){

  ldRanges <- make_ranges(ld$X1, ld$X2, ld$X3)
  ldRanges <- plyranges::mutate(ldRanges, locus=ld$X4)

  snpRanges <- make_ranges(seqname = cleaned.sumstats$chr,
                           start = cleaned.sumstats$pos,
                           end = cleaned.sumstats$pos)

  snpRanges <- plyranges::mutate(snpRanges, snp=cleaned.sumstats$snp)

  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
  snp.ld.block <- tibble::as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple ld-blocks due to edge of ld blocks
  cleaned.annot.sumstats <- dplyr::inner_join(cleaned.sumstats, snp.ld.block, 'snp')

  return(cleaned.annot.sumstats)
}


#' @title Assigns SNPs with annotations based on overlap
#' @param gwas a data frame or tibble of GWAS summary statistics
#' @param annotations annotation BED files
#'
#' @export
annotator <- function(gwas, annotations){

  snpRanges <- make_ranges(gwas$chr, gwas$pos, gwas$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=gwas$snp)

  for(f in annotations){

    name <- paste0(basename(f),'_d')
    curr <- rtracklayer::import(f, format='bed')
    subdf <- IRanges::subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)

    gwas <- dplyr::mutate(gwas, !!name := ifelse(snp %in% snpsIn,1,0))
  }
  return(gwas)
}

#' @title Annotations for causal SNPs (apply these after fine-mapping!)
#' @param gwas a data frame or tibble of GWAS summary statistics
#' @param annotations annotation BED files
#' @importFrom magrittr %>%
#' @export
annotator_merged <- function(gwas, annotations){

  snpRanges <- make_ranges(gwas$chr, gwas$pos, gwas$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=gwas$snp)
  gwas['annots'] <- ''

  for(f in annotations){

    curr <- rtracklayer::import(f, format='bed')
    subdf <- IRanges::subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)

    if(length(snpsIn)>0){
      curr <- gwas %>% pull(annots)
      curr <- curr[gwas$snp %in% snpsIn]
      delims <- rep(';', length(curr))
      delims[which(curr == '')] <- ''
      gwas[gwas$snp %in% snpsIn,"annots"] <- paste0(curr,delims,gsub(pattern = '.bed',replacement = '', x = basename(f)))
    }
  }
  return(gwas)
}

#' @title merge_bigsnp_gwas
#' @param gwas  GWAS summary statistics
#' @param bigSNP  bigSNP object
#' @importFrom magrittr %>%
#' @export
merge_bigsnp_gwas <- function(gwas, bigSNP){

  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')

  matched.gwas <- tibble::as_tibble(bigsnpr::snp_match(gwas,
                                                       snp_info,
                                                       strand_flip = TRUE,
                                                       match.min.prop = 0.1))

  matched.gwas <- matched.gwas %>%
    dplyr::rename(og_index = `_NUM_ID_.ss`) %>%
    dplyr::rename(bigSNP_index = `_NUM_ID_`) %>%
    dplyr::mutate(zscore = beta/se)

  return(matched.gwas)
}
