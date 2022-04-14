#' @title Clean GWAS summary statistics and adds metadata
#' @description Cleans GWAS summary statistics and adds metadata
#' @param sumstats_file file containing raw summary statistics.
#' Coordinates should be hg39/b37; the pipeline does not support hg38.
#' @param cols.to.keep character vector of the following columns:
#' chr, position, allele1, allele2, beta, se, unique id, pvalue
#' @param bigSNP bigSNP object from \code{bigsnpr}.
#' @param LD_Blocks Reference LD blocks
#' @return Cleaned summary statistics + LD block of every SNP,
#' as well as its index in the reference panel of genotypes
#' @export
run_gwas_cleaner <- function(sumstats_file, cols.to.keep, bigSNP, LD_Blocks){

  cat('Loading summary statistics...\n')
  if(is.character(sumstats_file) & length(sumstats_file) == 1){
    sumstats <- vroom::vroom(sumstats_file, col_names = TRUE, show_col_types = FALSE)
  }

  cat('Cleaning summary statistics...\n')
  cleaned_sumstats <- clean_sumstats(sumstats, cols.to.keep)

  cat('Matching to reference panel...\n')
  cleaned_sumstats <- merge_bigsnp_gwas(cleaned_sumstats, bigSNP = bigSNP)

  cat('Assining SNPs to LD blocks...\n')
  if(missing(LD_Blocks)){
    cat('No LD blocks supplied. Using the included 1KG European LD blocks.')
    data('Euro_LD_Chunks', package='Mapgen')
  }
  cleaned_sumstats <- assign_locus_snp(cleaned.sumstats = cleaned_sumstats,
                                       ld = LD_Blocks)

  return(cleaned_sumstats)

}

#' @title Clean summary statistics
#' @param sumstats a tibble or data frame containing raw summary statistics.
#' @param cols.to.keep columns to keep in the data frame
#'
#' @export
clean_sumstats <- function(sumstats, cols.to.keep){

  stopifnot(!is.null(sumstats))
  stopifnot(length(cols.to.keep) == 8)

  chr <- cols.to.keep[1]
  pos <- cols.to.keep[2]
  beta <- cols.to.keep[3]
  se <- cols.to.keep[4]
  a0 <- cols.to.keep[5]
  a1 <- cols.to.keep[6]
  rs <- cols.to.keep[7]
  pval <- cols.to.keep[8]

  # keep SNPs in 1kg
  # sumstats <- inner_join(sumstats, snps.to.keep, by=rs)

  # Extract relevant columns
  clean.sumstats <- sumstats[ ,c(chr, pos, beta, se, a0, a1, rs, pval)]
  colnames(clean.sumstats) <- c('chr','pos','beta','se','a0','a1','snp','pval')

  # drop XY chromosomes
  cat('Remove chrX and chrY...\n')
  clean.sumstats <- clean.sumstats[!(clean.sumstats$chr %in% c("X","Y")), ]
  # make chromosomes integers
  clean.sumstats$chr <- as.integer(clean.sumstats$chr)

  # Compute Zscores
  cat('Compute z-scores ...\n')
  zscore <- clean.sumstats$beta/clean.sumstats$se
  clean.sumstats['zscore'] <- zscore
  clean.sumstats <- clean.sumstats[!is.na(zscore),]

  # convert alleles to upper case
  clean.sumstats$a0 <- toupper(clean.sumstats$a0)
  clean.sumstats$a1 <- toupper(clean.sumstats$a1)

  # Keep SNPs only, no indels
  cat('Remove indels ...\n')
  nucs <- c('A','C','T','G')
  bola1 <- (clean.sumstats$a0 %in% nucs)
  bola2 <- (clean.sumstats$a1 %in% nucs)
  clean.sumstats <- clean.sumstats[bola1 & bola2,]

  # Sort by chromosome and position
  clean.sumstats <- clean.sumstats[order(clean.sumstats$chr, clean.sumstats$pos), ]

  # drop duplicate SNPs
  cat('Remove duplicate SNPs ...\n')
  chrpos <- paste0(clean.sumstats$chr, '_', clean.sumstats$pos)
  clean.sumstats <- clean.sumstats[!duplicated(chrpos), ]

  return(clean.sumstats)
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


#' @title Assign SNPs with annotations based on overlap
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
