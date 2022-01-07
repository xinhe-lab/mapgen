
#' @title Map SNPs to genes and assign weights
#'
#' @param finemap.gr Genomic Ranges of finemapping result.
#' @param genomic.annots A list of genomic annotations.
#' @param intron.mode If FALSE, do not assign intronic SNPs to genes containing the introns.
#' @param active_promoter_method Method to define active promoters.
#' @param enhancer_loop_method Enhancer loop method
#' @param c.dist A scaling number used for computing weight based on SNP-gene distance. Weight = exp(-dist/c). Default = 1e5 (100kb).
#' @param dist.to Compute distance to TSS, gene midpoint, or end (right point).
#' @param cols.to.keep columns to keep in the SNP gene weights
#' @import GenomicRanges
#' @import tidyverse
#'
#' @export
compute_gene_pip_by_hierarchies <- function(finemap.gr,
                                            genomic.annots,
                                            intron.mode = TRUE,
                                            active_promoter_method = c('OCRs', 'pcHiC', 'ABC', 'distance','none'),
                                            enhancer_loop_method = 'ABC.pcHiC.nearby20kb',
                                            c.dist = 1e5,
                                            dist.to = c('end', 'tss', 'midpoint'),
                                            cols.to.keep = c('snp','chr','pos','ref','alt','locus','pip','gene_name', 'category', 'weight')) {

  cat('Map SNPs to genes and assign weights ...\n')

  active_promoter_method <- match.arg(active_promoter_method)
  dist.to <- match.arg(dist.to)

  ## Define annotation categories

  # Exon and active promoters

  # Define active promoters
  if(active_promoter_method == 'OCRs' || is.null(genomic.annots$active.promoters)){
    cat('Define active promoter by overlapping with OCRs ...\n')
    genomic.annots$active.promoters <- IRanges::subsetByOverlaps(genomic.annots$promoters, genomic.annots$OCRs_hg19, minoverlap = 100)
  }else if(active_promoter_method == 'pcHiC'){
    cat('Define active promoter by promoters of pcHiC target ...\n')
    genomic.annots$active.promoters <- genomic.annots$promoters[genomic.annots$promoters$gene_name %in% genomic.annots$pcHiC$gene_name]
  }else if(active_promoter_method == 'ABC'){
    cat('Define active promoter by promoters of ABC target ...\n')
    genomic.annots$active.promoters <- genomic.annots$promoters[genomic.annots$promoters$gene_name %in% genomic.annots$ABC$gene_name]
  }else if(active_promoter_method == 'distance'){
    cat('Define active promoter by distance (2kb) ...\n')
    genomic.annots$active.promoters <- genomic.annots$promoters
  }else if(active_promoter_method == 'none'){
    cat('Do not use promoters...\n')
    genomic.annots$active.promoters <- NULL

  }

  if(active_promoter_method == 'none'){
    exons_active_promoters <- list(exons=genomic.annots$exons)
  }else{
    exons_active_promoters <- list(exons=genomic.annots$exons, active.promoters=genomic.annots$active.promoters)
  }

  # Enhancer loops
  if(grepl('base', enhancer_loop_method, ignore.case = T) || is.na(enhancer_loop_method)){
    cat('No enhancer loops. \n')
    enhancer_loops <- NULL
  }else{
    enhancer_loops <- list()
    if(grepl('ABC', enhancer_loop_method, ignore.case = T) && !is.null(genomic.annots$ABC)){
      cat('Include ABC scores in enhancer loops ... \n')
      enhancer_loops$ABC <- genomic.annots$ABC[, c('gene_name')]
    }

    if(grepl('pcHiC', enhancer_loop_method, ignore.case = T) && !is.null(genomic.annots$pcHiC)){
      cat('Include pcHiC in enhancer loops ... \n')
      enhancer_loops$pcHiC <- genomic.annots$pcHiC[, c('gene_name')]
    }

    if(grepl('coacc', enhancer_loop_method, ignore.case = T) && !is.null(genomic.annots$coacc)){
      cat('Include coacc in enhancer loops ... \n')
      enhancer_loops$coacc <- genomic.annots$coacc[, c('gene_name')]
    }

    if(grepl('nearby', enhancer_loop_method, ignore.case = T)){
      if(grepl('nearby20kb', enhancer_loop_method, ignore.case = T)){
        cat('Include enhancers with nearby promoters (20kb) in enhancer loops ... \n')
        enhancer_loops$nearby20kb <- genomic.annots$enhancer_nearby_promoter_20kb[, c('gene_name')]
      }else if(grepl('nearby10kb', enhancer_loop_method, ignore.case = T)){
        cat('Include enhancers with nearby promoters (10kb) in enhancer loops ... \n')
        enhancer_loops$nearby10kb <- genomic.annots$enhancer_nearby_promoter_10kb[, c('gene_name')]
      }
    }
  }

  # Introns and UTRs
  if(intron.mode) {
    cat('Assign intronic SNPs to genes containing the introns ...\n')
    intron_utrs <- list(introns=genomic.annots$introns, UTRs=genomic.annots$UTRs)
  }else{
    cat('Do not directly assign intronic SNPs to genes containing the introns ...\n')
    intron_utrs <- list(UTRs=genomic.annots$UTRs)
  }

  ## Overlap finemapped SNPs with each of our annotations
  cat('Overlapping finemapped SNPs with annotations ...\n')

  # Hierarchy level 1: first assign SNPs in exons and active promoters.
  cat('Hierarchy level 1: Assign SNPs in exons and active promoters ...\n')
  exons_active_promoters_overlap <- lapply(exons_active_promoters, function(x){plyranges::join_overlap_inner(x, finemap.gr)})
  snps.in <- unique(unlist(GenomicRanges::GRangesList(exons_active_promoters_overlap))$snp)
  finemap.gr <- finemap.gr[!(finemap.gr$snp %in% snps.in),]

  # Hierarchy level 2A: assign SNPs in enhancers to linked genes through enhancer loops.
  if(length(enhancer_loops) > 0){
    cat('Hierarchy level 2: Assign SNPs in enhancers to linked genes through enhancer loops ...\n')
    enhancer_loops_overlap <- lapply(enhancer_loops, function(x){plyranges::join_overlap_inner(x, finemap.gr)})
    snps.in <- unique(unlist(GenomicRanges::GRangesList(enhancer_loops_overlap))$snp)
    finemap.gr <- finemap.gr[!(finemap.gr$snp %in% snps.in),]
  }else{
    cat('Hierarchy level 2: Skipped. No enhancer loops ... \n')
    enhancer_loops_overlap <- NULL
  }

  # Hierarchy level 3: assign the other SNPs in enhancer regions to genes with distance weights.
  cat('Hierarchy level 3: Use enhancer mode. Assign SNPs in enhancer regions to genes with distance weights ...\n')
  finemap_in_enhancer_regions <- IRanges::subsetByOverlaps(finemap.gr, genomic.annots$enhancer_regions)
  enhancer_snps_genes_by_distance <- list(enhancer.regions = gene_by_distance(snp.gr = finemap_in_enhancer_regions, promoters.gr = genomic.annots$promoters, c.dist = c.dist, dist.to = dist.to))
  finemap.gr <- finemap.gr[!(finemap.gr$snp %in% finemap_in_enhancer_regions$snp),]

  # Hierarchy level 4: assign SNPs in introns/UTRs to genes containing the introns/UTRs.
  cat('Hierarchy level 4: Assign SNPs in introns/UTRs... \n')
  intron_utr_overlap <- lapply(intron_utrs, function(x){plyranges::join_overlap_inner(x, finemap.gr)})
  snps.in <- unique(unlist(GenomicRanges::GRangesList(intron_utr_overlap))$snp)
  finemap.gr <- finemap.gr[!(finemap.gr$snp %in% snps.in),]

  # Hierarchy level 5: assign the rest SNPs (intergenic) to genes with distance weights.
  cat('Hierarchy level 5: Assign SNPs (intergenic) to genes with distance weights ...\n')
  intergenic_snps_genes_by_distance <- list(intergenic = gene_by_distance(snp.gr = finemap.gr, promoters.gr = genomic.annots$promoters, c.dist = c.dist, dist.to = dist.to))
  snp.overlap.list <- c(exons_active_promoters_overlap, enhancer_loops_overlap, enhancer_snps_genes_by_distance, intron_utr_overlap, intergenic_snps_genes_by_distance)

  ## Assign weights
  cat('Assign weights to SNP-gene pairs ... \n')
  mat.list <- lapply(names(snp.overlap.list), function(x){
    snp.overlap.list[[x]] %>%
      as_tibble() %>%
      dplyr::mutate(category = x) %>%
      dplyr::distinct(gene_name, snp, category, .keep_all = T)
  })

  weights.mat <- Reduce(bind_rows, mat.list) %>% dplyr::select(all_of(cols.to.keep))
  weights.mat$weight <- ifelse( weights.mat$category %in% c(names(exons_active_promoters_overlap), names(enhancer_loops_overlap), names(intron_utr_overlap)),
                               1, weights.mat$weight )

  ## Combine results from different categories
  weights.mat <- weights.mat %>%
    dplyr::arrange(category) %>%
    dplyr::group_by(snp, gene_name) %>%
    dplyr::mutate(category = paste(category, collapse = ',')) %>%
    dplyr::distinct(snp, gene_name, .keep_all = TRUE)

  cat('SNP-gene pairs in each category: \n')
  print(table(weights.mat$category))

  cat('Compute gene PIP ... \n')
  # For each SNP, distribute PIP of a SNP to all linked genes. So weights of a SNP get normalized across genes (each SNP's weights of all genes should sum to 1)
  normalized.weights.mat <- weights.mat %>%
    dplyr::distinct(snp, gene_name, .keep_all = T) %>%
    dplyr::group_by(snp) %>%
    dplyr::mutate(frac_pip = weight/sum(weight))

  # Gene level:
  # For each gene, aggregate the fractional PIPs from all linked SNPs
  # gene PIP = sum of the fractional PIPs (pip * frac_pip) the gene received from all linked SNPs
  gene.pip.mat <- normalized.weights.mat %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(gene_pip = sum(pip * frac_pip)) %>%
    dplyr::arrange(gene_name)
  gene.pip.mat <- gene.pip.mat %>%
    dplyr::select(snp, gene_name, frac_pip, gene_pip) %>%
    dplyr::ungroup()

  snp.gene.pip.mat <- normalized.weights.mat %>% dplyr::left_join(., gene.pip.mat)
  snp.gene.pip.mat <- snp.gene.pip.mat[!is.na(snp.gene.pip.mat$gene_name),]

  return(snp.gene.pip.mat)
}


#' @title Extract gene level result from SNP level gene mapping result
#'
#' @param snp.gene.pip.res A data frame of SNP level gene mapping result
#' @param gene.annots a GRange object of gene annotations
#' @import tidyverse
#' @export
extract_gene_level_result <- function(snp.gene.pip.res, gene.annots) {
  cat('Extract gene level result ...\n')

  snp.gene.pip.res <- snp.gene.pip.res[!is.na(snp.gene.pip.res$gene_name),]

  gene.locations <- as.data.frame(gene.annots)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]
  gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))

  m <- match(snp.gene.pip.res$gene_name, gene.locations$gene_name)
  snp.gene.pip.res$gene_chr <- gene.locations[m, 'seqnames']
  snp.gene.pip.res$gene_pos <- gene.locations[m, 'tss']

  gene.pip.res <- snp.gene.pip.res[, c('gene_chr', 'gene_pos', 'gene_name', 'gene_pip')] %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::rename(chr = gene_chr, pos = gene_pos)

  return(gene.pip.res)

}

#' @title Get gene credible sets from SNP level gene mapping table
#'
#' @param snp.gene.pip.mat A data frame of SNP level gene mapping table
#' @param by.locus Logical, if TRUE, get credible sets based on locus level gene PIP
#' @param gene.cs.percent.thresh percentage threshold for gene credible sets
#' @import tidyverse
#' @export
gene_cs <- function(snp.gene.pip.mat,
                    by.locus = TRUE,
                    gene.cs.percent.thresh = 0.8){

  # add locus level gene PIP
  # For each locus - gene pair, sum over the fractional PIPs for SNPs in the locus and linked to the gene
  snp.locus.gene.pip.mat <- snp.gene.pip.mat %>%
    dplyr::group_by(locus, gene_name) %>%
    dplyr::mutate(locus_gene_pip = sum(pip * frac_pip)) %>% dplyr::ungroup()

  # simplify to get locus, gene_name, locus_gene_pip, and gene_pip
  locus.gene.pip.df <- snp.locus.gene.pip.mat %>%
    dplyr::select(locus, gene_name, gene_pip, locus_gene_pip) %>%
    dplyr::distinct(locus, gene_name, .keep_all=TRUE)

  # check if gene PIP is equal to the sum of gene-locus PIP
  for(gene in unique(locus.gene.pip.df$gene_name)){
    if(!all.equal(sum(locus.gene.pip.df$locus_gene_pip[locus.gene.pip.df$gene_name == gene]),
                  locus.gene.pip.df$gene_pip[locus.gene.pip.df$gene_name == gene][1])){
      cat('gene_pip of ', gene, ' is not equal to the sum of locus_gene_pip! \n')
    }
  }

  if(by.locus){
    # for each locus, keep the genes with gene locus PIP cumsum > 0.8
    gene.cumsum.df <- locus.gene.pip.df %>%
      dplyr::group_by(locus) %>%
      dplyr::arrange(desc(locus_gene_pip)) %>%
      dplyr::mutate(gene_pip_csum = cumsum(locus_gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.percent.thresh)[1])

    # create gene cs table
    gene.cs.df <- gene.cumsum.df %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(gene_cs = paste0(gene_name, collapse=','),
                gene_cs_locus_pip = paste(paste0(gene_name, '(',round(locus_gene_pip,3),')'), collapse=','),
                top_gene = gene_name[1],
                top_locus_gene_pip = locus_gene_pip[1],
                top_gene_pip = gene_pip[1])
  }else{
    # for each locus, keep the genes with gene locus PIP cumsum > 0.8
    gene.cumsum.df <- locus.gene.pip.df %>%
      dplyr::group_by(locus) %>%
      dplyr::arrange(desc(gene_pip)) %>%
      dplyr::mutate(gene_pip_csum = cumsum(gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.percent.thresh)[1])

    # create gene cs table
    gene.cs.df <- gene.cumsum.df %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(gene_cs = paste0(gene_name, collapse=','),
                gene_cs_pip = paste(paste0(gene_name, '(',round(gene_pip,3),')'), collapse=','),
                top_gene = gene_name[1],
                top_gene_pip = gene_pip[1])
  }

  return(list(gene.cumsum.df = gene.cumsum.df,
              gene.cs.df = gene.cs.df,
              locus.gene.pip.df = locus.gene.pip.df))
}



#' @title Find all genes (promoters) within 1MB genes and assign weights based on distance
#'
#' @param snp.gr Genomic Ranges with the SNP locations.
#' @param promoters.gr Genomic Ranges with the gene promoter locations.
#' @param c.dist A scaling number used for computing weight based on SNP-gene distance.
#' Weight = exp(-dist/c). Default = 1e5 (100kb).
#'
gene_by_distance <- function(snp.gr, promoters.gr, c.dist = 1e5,
                             dist.to = c('tss', 'midpoint', 'end')){

  # get all promoters within 1MB
  res <- plyranges::join_overlap_inner(x = promoters.gr,
                                       y = snp.gr,
                                       maxgap = 1e6)

  if(dist.to == 'tss'){
    res$gene_pos <- res$tss
  }else if(dist.to == 'midpoint'){
    res$gene_pos <- round((start(res) + end(res))/2)
  }else{
    res$gene_pos <- end(res)
  }

  res$distance <- abs(res$pos - res$gene_pos)
  res$weight <- exp(-res$distance / as.numeric(c.dist))

  return(res)
}


#' @title Find nearest genes (by distance to TSS or gene body)
#'
#' @param snp.gr GenomicRange object of SNP information
#' @param gene.gr GenomicRange object of gene information
#' @param dist.to Find nearest genes by distance to TSS or gene body
#' @import tidyverse
#' @export
#'
find_nearest_genes <- function(snp.gr, gene.gr,
                               dist.to = c('tss', 'genebody')){

  dist.to <- match.arg(dist.to)

  GenomeInfoDb::seqlevelsStyle(snp.gr) <- 'UCSC'
  GenomeInfoDb::seqlevelsStyle(gene.gr) <- 'UCSC'

  gene.locations <- as.data.frame(gene.gr)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]

  snp.gr$nearest_gene <- NA

  if(dist.to == 'tss'){
    gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))
    gene.locations$start <- gene.locations$end <- gene.locations$tss
    gene.locations.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations, keep.extra.columns = T)

    snp.nearest.gene.idx <- GenomicRanges::nearest(snp.gr, gene.locations.gr)
    snp.gr$nearest_gene <- gene.locations.gr$gene_name[snp.nearest.gene.idx]
  }else if(dist.to == 'genebody'){
    gene.locations.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations, keep.extra.columns = T)
    snp.nearest.gene.hits <- as.data.frame(nearest(snp.gr, gene.locations.gr, select = 'all'))
    colnames(snp.nearest.gene.hits) <- c('snp_idx', 'gene_idx')
    snp.nearest.gene.hits$snp <- snp.gr$snp[snp.nearest.gene.hits$snp_idx]
    snp.nearest.gene.hits$gene_name <- gene.locations.gr$gene_name[snp.nearest.gene.hits$gene_idx]

    snp.nearest.gene.df <- snp.nearest.gene.hits %>%
      dplyr::select(snp, gene_name) %>%
      dplyr::group_by(snp) %>%
      dplyr::summarise(nearestGene = paste(gene_name, collapse = ','))

    snp.nearest.gene.df$nearestGene <- as.character(snp.nearest.gene.df$nearestGene)

    snp.gr$nearest_gene <- snp.nearest.gene.df$nearestGene[match(snp.gr$snp, snp.nearest.gene.df$snp)]

  }

  return(snp.gr)
}


