
make_ranges <- function(seqname, start, end){
  return(GenomicRanges::GRanges(seqnames = seqname,
                                ranges = IRanges::IRanges(start = start, end = end)))
}


#' @title Assigns each SNP to one ld-block
#' @param cleaned.sumstats
#' @param ld
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @export
assign.locus.snp <- function(cleaned.sumstats, ld){

  ldRanges <- make_ranges(ld$X1, ld$X2, ld$X3)
  ldRanges <- plyranges::mutate(ldRanges, locus=ld$X4)

  snpRanges <- make_ranges(seqname = cleaned.sumstats$chr,
                           start = cleaned.sumstats$pos,
                           end = cleaned.sumstats$pos)

  snpRanges <- plyranges::mutate(snpRanges, snp=cleaned.sumstats$snp)

  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
  snp.ld.block <- as_tibble(snp.ld.overlap@elementMetadata)
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

#' @title merge.bigsnp.gwas
#' @param gwas  GWAS summary statistics
#' @param bigSNP  bigSNP object
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom bigsnpr snp_match
#' @export
merge.bigsnp.gwas <- function(gwas, bigSNP){

  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')

  matched.gwas <- as_tibble(bigsnpr::snp_match(gwas,
                                               snp_info,
                                               strand_flip = T,
                                               match.min.prop = 0.1))

  matched.gwas <- matched.gwas %>%
    dplyr::rename(og_index = `_NUM_ID_.ss`) %>%
    dplyr::rename(bigSNP_index = `_NUM_ID_`) %>%
    dplyr::mutate(zscore = beta/se)

  return(matched.gwas)
}

