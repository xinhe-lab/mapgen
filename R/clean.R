#' @title RunCleaner
#' @description Cleans GWAS summary statistics and adds metadata
#' @param sumstats a tibble or data frame containing raw summary statistics. Coordinates should be hg39/b37; the pipeline does not support hg38.
#' @param ColsToKeep character vector of the following columns: chr, position, allele1, allele2, beta, se, unique id, pvalue
#' @return Cleaned summary statistics + LD block of every SNP, as well as its index in the reference panel of genotypes
#' @export
RunCleaner <- function(sumstats, ColsToKeep, bigSNP){

  print('Loading summary statistics...')

  if(is.character(sumstats) & length(sumstats) == 1){

   sumstats <- vroom::vroom(sumstats, col_names = TRUE)

  }

  print('Cleaning summary statistics..')
  cleaned_sumstats <- clean_sumstats(sumstats, ColsToKeep)

  print('Matching to reference panel...')
  cleaned_sumstats <- merge.bigsnp.gwas(cleaned_sumstats, bigSNP = bigSNP)

  print('Assining SNPs to LD blocks...')
  data('Euro_LD_Chunks', package='Mapgen')
  cleaned_sumstats <- assign.locus.snp(cleaned.sumstats = cleaned_sumstats, ld = LD_Blocks)

  print('Complete.')

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
  #sumstats <- inner_join(sumstats, snps.to.keep, by=rs)
  # Extract relevant columns
  clean.sumstats <- sumstats[ ,c(chr, pos, beta, se, a0, a1, rs, pval)]
  colnames(clean.sumstats) <- c('chr','pos','beta','se','a0','a1','snp', 'pval')

  # drop XY chromosomes
  clean.sumstats <- clean.sumstats[!(clean.sumstats$chr %in% c("X","Y")), ]
  print('X,Y dropped')
  # make chromosomes integers
  clean.sumstats$chr <- as.integer(clean.sumstats$chr)

  # Compute Zscores
  zscore <- clean.sumstats$beta/clean.sumstats$se
  clean.sumstats['zscore'] <- zscore
  clean.sumstats <- clean.sumstats[!is.na(zscore),]
  print('zscore computed')

  # convert alleles to upper case
  clean.sumstats$a0 <- toupper(clean.sumstats$a0)
  clean.sumstats$a1 <- toupper(clean.sumstats$a1)

  # Keep SNPs only, no indels
  nucs <- c('A','C','T','G')
  bola1 <- (clean.sumstats$a0 %in% nucs)
  bola2 <- (clean.sumstats$a1 %in% nucs)
  clean.sumstats <- clean.sumstats[bola1 & bola2,]

  print('indels dropped')
  # sort by chromosome and position
  clean.sumstats <- clean.sumstats[order(clean.sumstats$chr, clean.sumstats$pos), ]

  # drop duplicate SNPs
  chrpos <- paste0(clean.sumstats$chr, '_', clean.sumstats$pos)
  clean.sumstats <- clean.sumstats[!duplicated(chrpos), ]
  print('duplicate snps removed')

  return(clean.sumstats)
}
