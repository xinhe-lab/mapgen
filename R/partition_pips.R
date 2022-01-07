

#' @title Sum PIPs for each annotation category
#'
#' @param finemap.locus.gr fine mapping result for a locus
#' @param annots.list list of functional annotations
#' @param unassigned name of unassigned category, default: 'Intergenic'.
#'
#' @export
#'
sum_pip_annots <- function(finemap.locus.gr, annots.list, unassigned = 'Intergenic'){

  sum.pip.categories <- matrix(NA, nrow = 2, ncol = length(annots.list)+1)
  colnames(sum.pip.categories) <- c(names(annots.list), unassigned)
  SNPsIn <- NULL
  for(i in 1:length(annots.list)){
    finemap.locus.gr <- finemap.locus.gr[!finemap.locus.gr$snp %in% SNPsIn, ]
    overlaps <- (GenomicRanges::countOverlaps(finemap.locus.gr, annots.list[[i]], ignore.strand = TRUE) > 0)
    finemap.locus.in.gr <- finemap.locus.gr[overlaps,]
    Sum.PIPs <- sum(finemap.locus.in.gr$pip)
    SNPs.in.annot <- unique(finemap.locus.in.gr$snp)
    N.SNPs <- length(SNPs.in.annot)
    sum.pip.categories[,i] <- c(Sum.PIPs, N.SNPs)
    SNPsIn <- c(SNPsIn, SNPs.in.annot)
  }

  finemap.locus.gr <- finemap.locus.gr[!finemap.locus.gr$snp %in% SNPsIn, ]

  Sum.PIPs <- sum(finemap.locus.gr$pip)
  N.SNPs <- length(unique(finemap.locus.gr$snp))
  sum.pip.categories[,unassigned] <- c(Sum.PIPs, N.SNPs)

  return(sum.pip.categories)
}


#' @title Sum PIPs for genomic regions
#'
#' @param finemap.locus.gr a GRanges object of fine mapping result for a locus.
#' @param regions.gr a GRanges object of genomic regions.
#' @param type Sum PIPs for SNPs inside or outside genomic regions.
#'
#' @export
#'
sum_pip_regions <- function(finemap.locus.gr, regions.gr, type = c('inside', 'outside')){
  type <- match.arg(type)
  overlaps <- (GenomicRanges::countOverlaps(finemap.locus.gr, regions.gr, ignore.strand = TRUE) > 0)

  if(type == 'inside'){
    finemap.locus.inside.regions.gr <- finemap.locus.gr[overlaps]
    Sum.PIPs = sum(finemap.locus.inside.regions.gr$pip)
    N.SNPs = length(finemap.locus.inside.regions.gr)
  }else if(type == 'outside'){
    finemap.locus.outside.regions.gr <- finemap.locus.gr[!overlaps]
    Sum.PIPs = sum(finemap.locus.outside.regions.gr$pip)
    N.SNPs = length(finemap.locus.outside.regions.gr)
  }
  return(c(Sum.PIPs = Sum.PIPs, N.SNPs = N.SNPs))
}


#' @title Partition PIPs into functional categories
#'
#' @param finemap.gr fine mapping result
#' @param annots functional annotations
#' @param unassigned name of unassigned category, default: 'Intergenic'.
#' @export
#'
partition_pip_annots <- function(finemap.gr, annots, unassigned = 'Intergenic'){
  locus.list <- sort(unique(finemap.gr$locus))
  # sum PIP matrix, rows are finemapped loci, columns are annotation categories
  sum.pip.mat <- matrix(NA, nrow = length(locus.list), ncol = length(annots))
  rownames(sum.pip.mat) <- locus.list
  for(i in 1:length(locus.list)){
    finemap.locus.gr <- finemap.gr[finemap.gr$locus == locus]
    sum.pip.categories <- as.numeric(sum_pip_annots(finemap.locus.gr, annots, unassigned))
    names(sum.pip.categories) <- paste(rep(c(names(annots), unassigned),each = 2), c("Sum.PIPs", "N.SNPs"))
    sum.pip.mat[i,] <- sum.pip.categories
  }
  sum.pip.df <- data.frame(locus = locus.list, sum.pip.mat, check.names = FALSE)
  return(sum.pip.df)
}

#' @title Partition PIPs into functional annotation regions
#'
#' @param finemap.gr fine mapping result
#' @param annots functional annotation regions
#' @param unassigned name of unassigned category, default: 'Intergenic'.
#' @export
#'
partition_pip_regions <- function(finemap.gr, annots, unassigned = 'Intergenic'){
  locus.list <- sort(unique(finemap.gr$locus))
  merged.annots <- unlist(as(annots, "GRangesList"))
  # sum PIP matrix, rows are finemapped loci, columns are annotation categories
  sum.pip.mat <- matrix(NA, nrow = length(locus.list), ncol = length(annots))
  rownames(sum.pip.mat) <- locus.list
  for(i in 1:length(locus.list)){
    finemap.locus.gr <- finemap.gr[finemap.gr$locus == locus]
    sum.pip.categories <- sapply(annots, function(x){sum_pip_regions(finemap.locus.gr, x, type = 'inside')})
    sum.pip.others <- sum_pip_regions(finemap.locus.gr, merged.annots, type = 'outside')
    sum.pip.categories <- as.numeric(cbind(sum.pip.categories, others = sum.pip.others))
    names(sum.pip.categories) <- paste(rep(c(names(annots),unassigned),each = 2), c("Sum.PIPs", "N.SNPs"))
    sum.pip.mat[i,] <- sum.pip.categories
  }
  sum.pip.df <- data.frame(locus = locus.list, sum.pip.mat, check.names = FALSE)
  return(sum.pip.df)
}

