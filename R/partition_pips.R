
#' @title Partition PIPs into disjoint functional annotation regions
#'
#' @param finemapstats.gr a GRanges object of fine mapping summary statistics
#' @param annots.list a list of GRanges objects of disjoint functional annotations.
#' Note: the annotation regions need to be disjoint!
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
    sum.pip.categories <- sapply(annots.list, function(x){sum_pip_regions(finemapstats.locus.gr, x)})
    sum.pip.mat[i,] <- as.numeric(sum.pip.categories)
  }

  sum.pips <- as.data.frame(sum.pip.mat[, grep('Sum.PIPs', colnames(sum.pip.mat), value = TRUE)])
  names(sum.pips) <- gsub(' Sum.PIPs','', names(sum.pips))

  n.snps <- as.data.frame(sum.pip.mat[, grep('N.SNPs', colnames(sum.pip.mat), value = TRUE)])
  names(n.snps) <- gsub(' N.SNPs','', names(n.snps))

  return(list(sum.pips = sum.pips, n.snps = n.snps))
}

#' Sum PIPs for genomic regions
#'
sum_pip_regions <- function(finemapstats.locus.gr, regions.gr){
  overlaps <- (GenomicRanges::countOverlaps(finemapstats.locus.gr, regions.gr, ignore.strand = TRUE) > 0)
  finemapstats.locus.in.gr <- finemapstats.locus.gr[overlaps]
  return(c(Sum.PIPs = sum(finemapstats.locus.in.gr$pip),
           N.SNPs = length(finemapstats.locus.in.gr)))
}


#' @title Partition PIPs into functional annotation categories.
#'
#' @param finemapstats.gr a GRanges object of fine mapping summary statistics
#' @param annots.list a list of GRanges objects of ordered functional annotation categories,
#' if a SNP is in multiple annotation categories,
#' it will be assigned to the first ordered category.
#' @return a list with a data frame with the sum of PIPs and
#' a data frame with the number of SNPs included for each annotation category.
#' @export
#'
partition_pip_annots <- function(finemapstats.gr, annots.list){
  locus.list <- sort(unique(finemapstats.gr$locus))

  sum.pip.mat <- matrix(NA, nrow = length(locus.list), ncol = (length(annots.list)+1)*2)
  rownames(sum.pip.mat) <- locus.list
  colnames(sum.pip.mat) <- paste(rep(c(names(annots.list), "others"),each = 2), c('Sum.PIPs', 'N.SNPs'))
  for(i in 1:length(locus.list)){
    finemapstats.locus.gr <- finemapstats.gr[finemapstats.gr$locus == locus.list[i]]
    sum.pip.categories <- sum_pip_annots(finemapstats.locus.gr, annots.list)
    sum.pip.mat[i,] <- as.numeric(sum.pip.categories)
  }

  sum.pips <- as.data.frame(sum.pip.mat[, grep('Sum.PIPs', colnames(sum.pip.mat), value = TRUE)])
  names(sum.pips) <- gsub(' Sum.PIPs','', names(sum.pips))

  n.snps <- as.data.frame(sum.pip.mat[, grep('N.SNPs', colnames(sum.pip.mat), value = TRUE)])
  names(n.snps) <- gsub(' N.SNPs','', names(n.snps))

  return(list(sum.pips = sum.pips, n.snps = n.snps))
}


#' Sum PIPs for each annotation category,
#' assign SNPs to annotations with the priority of the earlier category
sum_pip_annots <- function(finemapstats.locus.gr, annots.list){

  sum.pip.categories <- matrix(NA, nrow = 2, ncol = length(annots.list)+1)
  colnames(sum.pip.categories) <- c(names(annots.list), "others")

  counted.SNPs <- NULL
  for(i in 1:length(annots.list)){
    finemapstats.locus.gr <- finemapstats.locus.gr[!finemapstats.locus.gr$snp %in% counted.SNPs, ]
    overlaps <- (GenomicRanges::countOverlaps(finemapstats.locus.gr, annots.list[[i]], ignore.strand = TRUE) > 0)
    finemap.locus.in.gr <- finemapstats.locus.gr[overlaps]

    SNPs.in.annot <- unique(finemap.locus.in.gr$snp)
    Sum.PIPs <- sum(finemap.locus.in.gr$pip)
    N.SNPs <- length(SNPs.in.annot)
    sum.pip.categories[,i] <- c(Sum.PIPs, N.SNPs)
    counted.SNPs <- c(counted.SNPs, SNPs.in.annot)
  }

  finemapstats.locus.gr <- finemapstats.locus.gr[!finemapstats.locus.gr$snp %in% counted.SNPs, ]
  Sum.PIPs <- sum(finemapstats.locus.gr$pip)
  N.SNPs <- length(unique(finemapstats.locus.gr$snp))
  sum.pip.categories[,"others"] <- c(Sum.PIPs, N.SNPs)

  return(sum.pip.categories)
}

