
#' @title Sum PIPs for genomic regions
#'
#' @param finemapstats.locus.gr a GRanges object of fine mapping result for a locus.
#' @param regions.gr a GRanges object of genomic regions.
#' @param type Sum PIPs for SNPs inside or outside genomic regions.
#' @return a vector with the sum of PIPs and the number of SNPs included.
#' @export
#'
sum_pip_regions <- function(finemapstats.locus.gr, regions.gr, type = c('inside', 'outside')){
  type <- match.arg(type)
  overlaps <- (GenomicRanges::countOverlaps(finemapstats.locus.gr, regions.gr, ignore.strand = TRUE) > 0)

  if(type == 'inside'){
    finemapstats.locus.inside.regions.gr <- finemapstats.locus.gr[overlaps]
    Sum.PIPs = sum(finemapstats.locus.inside.regions.gr$pip)
    N.SNPs = length(finemapstats.locus.inside.regions.gr)
  }else if(type == 'outside'){
    finemapstats.locus.outside.regions.gr <- finemapstats.locus.gr[!overlaps]
    Sum.PIPs = sum(finemapstats.locus.outside.regions.gr$pip)
    N.SNPs = length(finemapstats.locus.outside.regions.gr)
  }
  return(c(Sum.PIPs = Sum.PIPs, N.SNPs = N.SNPs))
}


#' @title Partition PIPs into functional annotation regions
#'
#' @param finemapstats.gr a GRanges object of fine mapping summary statistics
#' @param annots.list a list of GRanges objects of functional annotations
#' @return a list with a data frame with the sum of PIPs and
#' a data frame with the number of SNPs included for each annotation category.
#' @export
#'
partition_pip_regions <- function(finemapstats.gr, annots.list){
  locus.list <- sort(unique(finemapstats.gr$locus))

  sum.pip.mat <- matrix(NA, nrow = length(locus.list), ncol = length(annots.list)*2)
  rownames(sum.pip.mat) <- locus.list
  colnames(sum.pip.mat) <- paste(rep(names(annots.list),each = 2), c('Sum.PIPs', 'N.SNPs'))
  for(i in 1:length(locus.list)){
    finemapstats.locus.gr <- finemapstats.gr[finemapstats.gr$locus == locus.list[i]]
    sum.pip.categories <- sapply(annots.list, function(x){sum_pip_regions(finemapstats.locus.gr, x, type = 'inside')})
    sum.pip.mat[i,] <- as.numeric(sum.pip.categories)
  }

  sum.pips <- as.data.frame(sum.pip.mat[, grep('Sum.PIPs', colnames(sum.pip.mat), value = TRUE)])
  names(sum.pips) <- gsub(' Sum.PIPs','', names(sum.pips))

  n.snps <- as.data.frame(sum.pip.mat[, grep('N.SNPs', colnames(sum.pip.mat), value = TRUE)])
  names(n.snps) <- gsub(' N.SNPs','', names(n.snps))

  return(list(sum.pips = sum.pips, n.snps = n.snps))
}

