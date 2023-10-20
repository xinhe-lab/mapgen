
#' @title Compute gene PIPs based on fine-mapping result and functional annotations.
#'
#' @param finemapstats A GRanges object of fine-mapping result.
#' @param genomic.annots A list of GRanges objects of genomic annotations.
#' @param intron.mode Logical. If TRUE, assign intronic SNPs to genes containing the introns.
#' @param d0 A scaling parameter used for computing weight based on SNP-gene distance.
#' Weight = exp(-dist/d0). Default = 50000 (50kb).
#' @param exon.weight Weight for exon (and promoter) categories (default = 1).
#' @param loop.weight Weight for chromatin loop categories (default = 1).
#' @param cols.to.keep columns to keep in the SNP gene weights
#' @return A data frame of gene mapping result
#' @export
compute_gene_pip <- function(finemapstats,
                             genomic.annots,
                             intron.mode = FALSE,
                             d0 = 50000,
                             exon.weight = 1,
                             loop.weight = 1,
                             cols.to.keep = c(names(mcols(finemapstats)), 'gene_name', 'category', 'weight', 'frac_pip', 'gene_pip')) {

  cat('Map SNPs to genes and assign weights ...\n')

  # Define gene mapping categories

  ## Exons and promoters category
  exons_promoters <- genomic.annots$exons_promoters
  if(is.null(exons_promoters)){
    exons_promoters <- list(exons=genomic.annots$exons,
                            promoters=genomic.annots$promoters)
  }

  ## Enhancer loops category
  enhancer_loops <- genomic.annots$enhancer_loops

  ## Enhancer regions category
  enhancer_regions <- genomic.annots$enhancer_regions

  ## Introns and UTRs category
  if(intron.mode){
    introns_UTRs <- list(introns=genomic.annots$introns,
                         UTRs=genomic.annots$UTRs)
  }else{
    introns_UTRs <- list(UTRs=genomic.annots$UTRs)
  }

  # Assign SNPs to genes by hierarchical model

  ## Hierarchy level 1: first assign SNPs in exons and active promoters.
  if(!is.null(exons_promoters)){
    cat('Assign SNPs in exons and promoters ...\n')
    exons_promoters <- lapply(exons_promoters, function(x){x <- x[,'gene_name']})
    exons_promoters_assignment <- lapply(exons_promoters,
                                         function(x){plyranges::join_overlap_inner(x, finemapstats)})
    snps.in <- unique(unlist(GenomicRanges::GRangesList(exons_promoters_assignment))$snp)
    finemapstats <- finemapstats[!(finemapstats$snp %in% snps.in),]
  }else{
    cat('Category of exons/promoters not available! \n')
    exons_promoters_assignment <- NULL
  }

  ## Hierarchy level 2: assign SNPs in enhancers
  ## Hierarchy level 2A: assign SNPs in enhancers to linked genes through enhancer loops.
  if(!is.null(enhancer_loops) > 0){
    cat('Assign SNPs to linked genes through enhancer loops ...\n')
    cat('Enhancer loops include:', paste(names(enhancer_loops), collapse=', '), '\n')
    enhancer_loops <- lapply(enhancer_loops, function(x){x <- x[,'gene_name']})
    enhancer_loops_assignment <- lapply(enhancer_loops, function(x){plyranges::join_overlap_inner(x, finemapstats)})
    snps.in <- unique(unlist(GenomicRanges::GRangesList(enhancer_loops_assignment))$snp)
    finemapstats <- finemapstats[!(finemapstats$snp %in% snps.in),]
  }else{
    cat('Category of enhancer loops not available! \n')
    enhancer_loops_assignment <- NULL
  }

  ## Hierarchy level 2B: assign the rest of SNPs in enhancer regions to genes by distance weighting.
  if(!is.null(enhancer_regions)){
    cat('Assign SNPs in enhancer regions to genes by distance weighting ...\n')
    finemapstats_in_enhancer_regions <- IRanges::subsetByOverlaps(finemapstats, enhancer_regions)
    enhancer_regions_by_distance <- list(
      enhancer_regions = compute_distance_weight(finemapstats_in_enhancer_regions,
                                                 genomic.annots$promoters,
                                                 d0 = d0, type = 'promoter'))
    snps.in <- finemapstats_in_enhancer_regions$snp
    finemapstats <- finemapstats[!(finemapstats$snp %in% snps.in),]
  }else{
    cat('Category of enhancer regions not available! \n')
    enhancer_regions_by_distance <- NULL
  }

  ## Hierarchy level 3: assign SNPs in introns/UTRs (excluding enhancer regions).
  if(!is.null(introns_UTRs)){
    cat('Assign SNPs in UTRs (excluding enhancer regions) to the UTR genes ... \n')
    if(intron.mode){
      cat('Assign SNPs in introns (excluding enhancer regions) to genes containing the introns ...\n')
    }
    introns_UTRs <- lapply(introns_UTRs, function(x){x <- x[,'gene_name']})
    introns_UTRs_assignment <- lapply(introns_UTRs, function(x){plyranges::join_overlap_inner(x, finemapstats)})
    snps.in <- unique(unlist(GenomicRanges::GRangesList(introns_UTRs_assignment))$snp)
    finemapstats <- finemapstats[!(finemapstats$snp %in% snps.in),]
  }else{
    cat('Category of introns/UTRs not available! \n')
    introns_UTRs_assignment <- NULL
  }

  ## Hierarchy level 4: assign the rest of SNPs (intergenic) to genes by distance weighting.
  if(length(finemapstats) > 0){
    cat('Assign the rest of SNPs (intergenic) to genes by distance weighting ...\n')
    intergenic_by_distance <- list(
      intergenic = compute_distance_weight(finemapstats, genomic.annots$promoters,
                                           d0 = d0, type = 'promoter'))
  }else{
    cat('Category of intergenic SNPs not available! \n')
    intergenic_snps_genes_by_distance <- NULL
  }

  cat('Compute gene PIP ... \n')
  # Combine results from different categories
  snp.assignment.list <- c(exons_promoters_assignment,
                           enhancer_loops_assignment,
                           enhancer_regions_by_distance,
                           introns_UTRs_assignment,
                           intergenic_by_distance)

  # Assign weights to SNP-gene pairs
  mat.list <- lapply(names(snp.assignment.list), function(x){
    snp.assignment.list[[x]] %>% tibble::as_tibble() %>%
      dplyr::mutate(category = x) %>%
      dplyr::distinct(gene_name, snp, category, .keep_all = T)
  })

  weights.mat <- Reduce(dplyr::bind_rows, mat.list)

  # Assign weights for different categories.
  gene_categories <- c(names(exons_promoters_assignment),
                       names(introns_UTRs_assignment))

  loop_categories <- names(enhancer_loops_assignment)

  distance_categories <- c(names(enhancer_regions_by_distance),
                           names(intergenic_by_distance))

  # cat('weight =', exon.weight, 'for:', gene_categories, '\n')
  weights.mat <- weights.mat %>%
    dplyr::mutate(weight = ifelse(category %in% gene_categories, exon.weight, weight))

  # cat('weight =', loop.weight, 'for:', loop_categories, '\n')
  weights.mat <- weights.mat %>%
    dplyr::mutate(weight = ifelse(category %in% loop_categories, loop.weight, weight))

  # cat('distance weight for:', distance_categories, '\n')
  weights.mat <- weights.mat %>%
    dplyr::mutate(category = ifelse(category %in% distance_categories, 'distance', category)) %>%
    dplyr::arrange(category) %>%
    dplyr::group_by(snp, gene_name) %>%
    dplyr::mutate(category = paste(category, collapse = ',')) %>%
    dplyr::distinct(snp, gene_name, .keep_all = TRUE)

  # For each SNP, distribute PIP of a SNP to all linked genes.
  # So weights of a SNP get normalized across genes (each SNP's weights of all genes should sum to 1)
  normalized.weights.mat <- weights.mat %>%
    dplyr::distinct(snp, gene_name, .keep_all = TRUE) %>%
    dplyr::group_by(snp) %>%
    dplyr::mutate(frac_pip = weight/sum(weight))

  # For each gene, aggregate the fractional PIPs from all linked SNPs
  # gene PIP = sum(pip * frac_pip)
  gene.pip.mat <- normalized.weights.mat %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(gene_pip = sum(pip * frac_pip)) %>%
    dplyr::arrange(gene_name)

  gene.pip.mat <- gene.pip.mat %>%
    dplyr::select(snp, gene_name, frac_pip, gene_pip) %>%
    dplyr::ungroup()

  gene.mapping.res <- suppressMessages(normalized.weights.mat %>% dplyr::left_join(., gene.pip.mat))
  gene.mapping.res <- gene.mapping.res[!is.na(gene.mapping.res$gene_name),]

  gene.mapping.res <- gene.mapping.res %>%
    dplyr::select(all_of(cols.to.keep)) %>%
    as.data.frame()

  cat("Done.")

  return(gene.mapping.res)
}


#' @title Extract gene-level result from SNP-level gene mapping result
#'
#' @param gene.mapping.res A data frame of SNP-level gene mapping result
#' @param gene.annots a GRanges object of gene annotations
#' @return A data frame of gene-level view of gene mapping result
#' @export
extract_gene_level_result <- function(gene.mapping.res, gene.annots) {

  gene.mapping.res <- gene.mapping.res[!is.na(gene.mapping.res$gene_name),]
  genes_not_included <- setdiff(gene.mapping.res$gene_name, gene.annots$gene_name)
  if(length(genes_not_included) > 0){
    message(sprintf('Remove %d genes not included in gene.annots.\n', length(genes_not_included)))
    gene.mapping.res <- gene.mapping.res %>% dplyr::filter(!gene_name %in% genes_not_included)
  }

  gene.locations <- as.data.frame(gene.annots)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]
  gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))

  m <- match(gene.mapping.res$gene_name, gene.locations$gene_name)
  gene.mapping.res$gene_chr <- gene.locations[m, 'seqnames']
  gene.mapping.res$gene_pos <- gene.locations[m, 'tss']

  gene.pip.res <- gene.mapping.res[, c('gene_chr', 'gene_pos', 'gene_name', 'gene_pip')] %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::rename(chr = gene_chr, pos = gene_pos) %>%
    as.data.frame()

  return(gene.pip.res)

}

#' @title Obtain credible gene sets from SNP-level gene mapping table
#'
#' @param gene.mapping.res A data frame of SNP-level gene mapping table
#' @param by.locus Logical, if TRUE, get credible gene sets based on locus-level gene PIP,
#' If FALSE, get credible gene sets based on gene PIP.
#' @param gene.cs.coverage A number between 0 and 1 specifying desired coverage of each credible gene set
#' @return a data frame of credible gene set result.
#' Columns:
#' gene_cs: credible gene sets,
#' gene_cs_locus_pip: credible gene sets and corresponding locus-level gene PIPs.
#' gene_cs_pip: credible gene sets and corresponding gene PIPs.
#' top_gene:  genes with highest gene PIP at each locus.
#' top_locus_gene_pip: locus-level gene PIP for the top gene.
#' Locus-level gene PIP only includes SNPs within a locus, so this value may be lower than the gene PIP.
#' top_gene_pip:  gene PIP of the top gene.
#' @export
gene_cs <- function(gene.mapping.res,
                    by.locus = TRUE,
                    gene.cs.coverage = 0.8){

  # Get locus level gene PIP
  locus.gene.pip.df <- get_locus_level_gene_pip(gene.mapping.res)

  # check to make sure gene PIP is equal to the sum of gene-locus PIP
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
      dplyr::filter(sum(locus_gene_pip) >= gene.cs.coverage) %>%
      dplyr::arrange(desc(locus_gene_pip)) %>%
      dplyr::mutate(gene_pip_csum = cumsum(locus_gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.coverage)[1])

    gene.cs.df <- gene.cumsum.df %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(gene_cs = paste0(gene_name, collapse=','),
                       gene_cs_locus_pip = paste(paste0(gene_name, '(',round(locus_gene_pip,3),')'), collapse=','),
                       top_gene = gene_name[1],
                       top_locus_gene_pip = round(locus_gene_pip[1],3),
                       top_gene_pip = round(gene_pip[1],3))
  }else{
    # for each locus, keep the genes with gene PIP cumsum > 0.8
    gene.cumsum.df <- locus.gene.pip.df %>%
      dplyr::group_by(locus) %>%
      dplyr::filter(sum(gene_pip) >= gene.cs.coverage) %>%
      dplyr::arrange(desc(gene_pip)) %>%
      dplyr::mutate(gene_pip_csum = cumsum(gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.coverage)[1])

    gene.cs.df <- gene.cumsum.df %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(gene_cs = paste0(gene_name, collapse=','),
                       gene_cs_pip = paste(paste0(gene_name, '(',round(gene_pip,3),')'), collapse=','),
                       top_gene = gene_name[1],
                       top_gene_pip = round(gene_pip[1],3))
  }

  gene.cs.df <- as.data.frame(gene.cs.df)
  return(gene.cs.df)
}

#' @title Gene view summary table
#'
#' @param gene.mapping.res A data frame of gene mapping result
#' @param gene.pip.thresh Filter genes with gene PIP cutoff (default: 0.1)
#' @return A data frame of gene view summary of gene mapping result
#' @export
#'
gene_view_summary <- function(gene.mapping.res, gene.pip.thresh = 0.1){
  gene.view.df <- gene.mapping.res %>%
    dplyr::mutate(fractional_PIP = pip * frac_pip) %>%
    dplyr::select(gene_name, gene_pip, fractional_PIP) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(gene_pip = round(gene_pip[1], 3),
                     n_snps_frac_pip_2percent = sum(fractional_PIP > 0.02)) %>%
    dplyr::filter(gene_pip >= gene.pip.thresh)
  return(gene.view.df)
}

#' @title SNP view summary table
#'
#' @param gene.mapping.res A data frame of gene mapping result
#' @param gene.annots A data frame of gene annotations
#' @param finemapstats A GRange object of fine-mapping summary statistics
#' @param fractional.pip.thresh Filter SNPs with fractional PIP cutoff (default: 0.02)
#' @return A data frame of SNP view summary of gene mapping result
#' @export
#'
snp_view_summary <- function(gene.mapping.res, gene.annots, finemapstats, fractional.pip.thresh = 0.02){
  high.conf.snp.df <- gene.mapping.res %>% dplyr::filter(fractional_PIP >= fractional.pip.thresh)

  snp.gene <- high.conf.snp.df %>% dplyr::select(snp, pos, gene_name)

  gene.locs.df <- gene.annots %>% tibble::as_tibble()
  gene.locs.df$TSS <- ifelse(gene.locs.df$strand=='+', gene.locs.df$start, gene.locs.df$end)

  snp.gene.dist <- snp.gene %>%
    dplyr::left_join(., gene.locs.df, on = 'gene_name') %>%
    dplyr::mutate(dist = abs(TSS - pos)) %>%
    dplyr::select(snp, gene_name, dist)

  high.conf.snp.df <- dplyr::inner_join(high.conf.snp.df, snp.gene.dist, on = c('snp','gene_name'))

  # Add nearest gene (distance to gene body)
  GenomeInfoDb::seqlevelsStyle(finemapstats) <- 'UCSC'
  snp.nearest.gene.gr <- find_nearest_genes(finemapstats, gene.annots, dist.to = 'genebody')
  snp.nearest.gene.df <- as.data.frame(snp.nearest.gene.gr)[, c('snp', 'nearest_gene')]
  high.conf.snp.df <- high.conf.snp.df %>% dplyr::left_join(., snp.nearest.gene.df, on='snp')

  # SNP view table
  high.conf.snp.df <- high.conf.snp.df %>%
    dplyr::select(-weight, -frac_pip) %>%
    dplyr::rename(SNP = snp,
                  PIP = pip,
                  `Gene Linked` = gene_name,
                  `Gene PIP`=gene_pip,
                  `Link Method`=category,
                  `Distance to Gene`=dist,
                  `Nearest Gene` = nearest_gene) %>%
    dplyr::mutate(PIP = round(PIP, 3), `Gene PIP` = round(`Gene PIP`, 3)) %>%
    dplyr::arrange(chr, pos)

  return(high.conf.snp.df)
}

#' @title LD block view summary table
#'
#' @param gene.mapping.res A data frame of gene mapping result
#' @param finemapstats A GRange object of fine-mapping summary statistics
#' @param gene.cs.coverage A number between 0 and 1 specifying desired coverage of each credible gene set
#' @return A data frame of LD block view summary of gene mapping result
#' @export
#'
block_view_summary <- function(gene.mapping.res, finemapstats, gene.cs.coverage = 0.8){

  # Gene CS based on locus level gene PIP
  gene.cs.l <- gene_cs(gene.mapping.res, by.locus = TRUE, gene.cs.coverage = gene.cs.coverage)
  gene.cs.df <- gene.cs.l$gene.cs.df
  gene.cumsum.df <- gene.cs.l$gene.cumsum.df
  locus.gene.pip.df <- gene.cs.l$locus.gene.pip.df

  # Add nearest genes to LD blocks
  GenomeInfoDb::seqlevelsStyle(finemapstats) <- 'UCSC'
  finemapstats <- finemapstats[order(abs(finemapstats$zscore), decreasing = T), ]
  top.snps <- finemapstats[!duplicated(finemapstats$locus), ]
  nearest_genebody_genes <- find_nearest_genes(top.snps, gene.annots, dist.to = 'genebody')

  locus_topsnp_nearest_genes.df <- nearest_genebody_genes %>% tibble::as_tibble() %>%
    dplyr::select(locus, nearest_gene)

  block.view.df <- dplyr::left_join(gene.cs.df, locus_topsnp_nearest_genes.df, by = 'locus') %>%
    dplyr::left_join(., ldblocks.truegenes.df, by = 'locus') %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(`80% Gene Credible Set` =  paste0(unique(gene_cs), collapse=','),
                     `Top Genes` = paste0(unique(top_gene), collapse=','),
                     `Top Gene PIP (locus level)` = paste0(unique(round(top_locus_gene_pip,3)), collapse=','),
                     `Top Gene PIP` = paste0(unique(round(top_gene_pip,3)), collapse=','),
                     `Nearest Gene` = paste0(unique(nearest_gene), collapse = ','))

  return(block.view.df)
}

#' @title Get locus level gene PIP
#'
#' @param gene.mapping.res A data frame of SNP-level gene mapping result
#' @return A data frame of locus-level gene PIP result
#' @export
#'
get_locus_level_gene_pip <- function(gene.mapping.res){

  # For each locus - gene pair, sum over the fractional PIPs for SNPs
  # in the locus and linked to the gene
  snp.locus.gene.pip.mat <- gene.mapping.res %>%
    dplyr::group_by(locus, gene_name) %>%
    dplyr::mutate(locus_gene_pip = sum(pip * frac_pip)) %>%
    dplyr::ungroup()

  # simplify to get locus, gene_name, locus_gene_pip, and gene_pip
  locus.gene.pip.df <- snp.locus.gene.pip.mat %>%
    dplyr::select(locus, gene_name, gene_pip, locus_gene_pip) %>%
    dplyr::distinct(locus, gene_name, .keep_all=TRUE) %>%
    as.data.frame()

  return(locus.gene.pip.df)
}

#' @title Find the nearest genes for top SNPs in each locus.
#'
#' @param top.snps A GRanges object of the GWAS summary statistics for the top SNPs
#' @param genes A GRanges object of gene information
#' @param dist.to Find nearest genes by distance to gene body or TSS
#' @param cols.to.keep columns to keep in the result
#' @return A data frame with SNP locations and nearest genes.
#' @export
#'
find_nearest_genes <- function(top.snps,
                               genes,
                               dist.to = c('genebody', 'tss'),
                               cols.to.keep = c('snp','chr','pos','nearest_gene')){

  dist.to <- match.arg(dist.to)

  GenomeInfoDb::seqlevelsStyle(top.snps) <- 'UCSC'
  GenomeInfoDb::seqlevelsStyle(genes) <- 'UCSC'

  gene.locations <- as.data.frame(genes)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]

  top.snps$nearest_gene <- NA

  if(dist.to == 'tss'){
    gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))
    gene.locations.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations,
                                                                 start.field = 'tss',
                                                                 end.field = 'tss',
                                                                 keep.extra.columns = T)
    snp.nearest.gene.idx <- GenomicRanges::nearest(top.snps, gene.locations.gr)
    top.snps$nearest_gene <- gene.locations.gr$gene_name[snp.nearest.gene.idx]
  }else if(dist.to == 'genebody'){
    gene.locations.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations, keep.extra.columns = T)
    snp.nearest.gene.hits <- as.data.frame(nearest(top.snps, gene.locations.gr, select = 'all'))
    colnames(snp.nearest.gene.hits) <- c('snp_idx', 'gene_idx')
    snp.nearest.gene.hits$snp <- top.snps$snp[snp.nearest.gene.hits$snp_idx]
    snp.nearest.gene.hits$gene_name <- gene.locations.gr$gene_name[snp.nearest.gene.hits$gene_idx]

    snp.nearest.gene.df <- snp.nearest.gene.hits %>%
      dplyr::select(snp, gene_name) %>%
      dplyr::group_by(snp) %>%
      dplyr::summarise(nearest_gene = paste(gene_name, collapse = ','))

    snp.nearest.gene.df$nearest_gene <- as.character(snp.nearest.gene.df$nearest_gene)
    top.snps$nearest_gene <- snp.nearest.gene.df$nearest_gene[match(top.snps$snp, snp.nearest.gene.df$snp)]
  }

  return(as.data.frame(top.snps)[,cols.to.keep])
}

# Find genes around SNPs and assign weights based on distance
compute_distance_weight <- function(snp.ranges, gene.ranges,
                                    d0 = 50000, maxgap = 1e6,
                                    type = c('promoter', 'gene')){

  if(is.null(gene.ranges$tss)){
    gene.ranges$tss <- get_tss(gene.ranges, type = type)
  }

  snp.ranges$pos <- GenomicRanges::start(snp.ranges)

  res <- plyranges::join_overlap_inner(gene.ranges, snp.ranges, maxgap = maxgap)
  res$distance <- abs(res$pos - res$tss)
  res$weight <- exp(-res$distance / as.numeric(d0))

  return(res)
}


