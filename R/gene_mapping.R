
#' @title Compute gene PIPs based on fine-mapping result and functional annotations.
#'
#' @param finemap.gr a GRanges object of fine-mapping result.
#' @param genomic.annots A list of GRanges objects of genomic annotations.
#' @param enhancer.loop.method Enhancer loop method
#' @param intron.mode Logical. If TRUE, assign intronic SNPs to genes containing the introns.
#' @param c.dist A scaling number used for computing weight based on SNP-gene distance. Weight = exp(-dist/c). Default = 50000 (50kb).
#' @param cols.to.keep columns to keep in the SNP gene weights
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @import GenomicRanges
#' @return A data frame of SNP-level view of gene mapping result
#' @export
compute_gene_pip <- function(finemap.gr,
                             genomic.annots,
                             enhancer.loop.method = 'ABC.pcHiC.nearby20kb',
                             intron.mode = FALSE,
                             c.dist = 50000,
                             cols.to.keep = c(names(mcols(finemap.gr)), 'gene_name', 'category', 'weight', 'frac_pip', 'gene_pip')) {

  cat('Map SNPs to genes and assign weights ...\n')

  # Define annotation categories

  ## Exon and active promoters category
  if(!is.null(genomic.annots$active_promoters)){
    exons_active_promoters <- list(exons=genomic.annots$exons, active_promoters=genomic.annots$active_promoters)
  }else{
    exons_active_promoters <- list(exons=genomic.annots$exons)
  }

  ## Enhancer loops category
  enhancer_loops <- list()
  if(grepl('ABC', enhancer.loop.method, ignore.case = T) && !is.null(genomic.annots$ABC)){
    cat('Include ABC scores in enhancer loops ... \n')
    enhancer_loops$ABC <- genomic.annots$ABC[, c('gene_name')]
  }
  if(grepl('pcHiC', enhancer.loop.method, ignore.case = T) && !is.null(genomic.annots$pcHiC)){
    cat('Include pcHiC in enhancer loops ... \n')
    enhancer_loops$pcHiC <- genomic.annots$pcHiC[, c('gene_name')]
  }
  if(grepl('coacc', enhancer.loop.method, ignore.case = T) && !is.null(genomic.annots$coacc)){
    cat('Include coacc in enhancer loops ... \n')
    enhancer_loops$coacc <- genomic.annots$coacc[, c('gene_name')]
  }
  if(grepl('nearby', enhancer.loop.method, ignore.case = T)){
    if(grepl('nearby20kb', enhancer.loop.method, ignore.case = T)){
      cat('Include enhancers with nearby promoters (20kb) in enhancer loops ... \n')
      enhancer_loops$nearby20kb <- genomic.annots$enhancer_nearby_promoter_20kb[, c('gene_name')]
    }else if(grepl('nearby10kb', enhancer.loop.method, ignore.case = T)){
      cat('Include enhancers with nearby promoters (10kb) in enhancer loops ... \n')
      enhancer_loops$nearby10kb <- genomic.annots$enhancer_nearby_promoter_10kb[, c('gene_name')]
    }
  }

  ## Introns and UTRs category
  if(intron.mode) {
    intron_utrs <- list(introns=genomic.annots$introns, UTRs=genomic.annots$UTRs)
  }else{
    intron_utrs <- list(UTRs=genomic.annots$UTRs)
  }

  # Add TSS to promoters
  if(is.null(genomic.annots$promoters$tss)){
    plus_strand <- which(as.character(strand(genomic.annots$promoters)) == '+')
    minus_strand <- which(as.character(strand(genomic.annots$promoters)) == '-')
    genomic.annots$promoters$tss <- NA
    genomic.annots$promoters$tss[plus_strand] <- end(genomic.annots$promoters)[plus_strand]
    genomic.annots$promoters$tss[minus_strand] <- start(genomic.annots$promoters)[minus_strand]
  }

  # Assign SNPs to genes
  ## Hierarchy level 1: first assign SNPs in exons and active promoters.
  if(!is.null(exons_active_promoters)){
    cat('Assign SNPs in exons and active promoters ...\n')
    exons_active_promoters_overlap <- lapply(exons_active_promoters, function(x){plyranges::join_overlap_inner(x, finemap.gr)})
    snps.in <- unique(unlist(GenomicRanges::GRangesList(exons_active_promoters_overlap))$snp)
    finemap.gr <- finemap.gr[!(finemap.gr$snp %in% snps.in),]
  }else{
    exons_active_promoters_overlap <- NULL
  }

  ## Hierarchy level 2A: assign SNPs in enhancers to linked genes through enhancer loops.
  if(length(enhancer_loops) > 0){
    cat('Assign SNPs in enhancers to linked genes through enhancer loops ...\n')
    enhancer_loops_overlap <- lapply(enhancer_loops, function(x){plyranges::join_overlap_inner(x, finemap.gr)})
    snps.in <- unique(unlist(GenomicRanges::GRangesList(enhancer_loops_overlap))$snp)
    finemap.gr <- finemap.gr[!(finemap.gr$snp %in% snps.in),]
  }else{
    enhancer_loops_overlap <- NULL
  }

  ## Hierarchy level 2B: assign SNPs in enhancer regions to genes by distance weighting.
  if(!is.null(genomic.annots$enhancer_regions)){
    cat('Assign SNPs in enhancer regions to genes by distance weighting ...\n')
    finemap_in_enhancer_regions <- IRanges::subsetByOverlaps(finemap.gr, genomic.annots$enhancer_regions)
    enhancer_snps_genes_by_distance <- list(enhancer.regions = gene_by_distance(snps.gr = finemap_in_enhancer_regions,
                                                                                promoters.gr = genomic.annots$promoters,
                                                                                c.dist = c.dist))
    finemap.gr <- finemap.gr[!(finemap.gr$snp %in% finemap_in_enhancer_regions$snp),]
  }else{
    enhancer_snps_genes_by_distance <- NULL
  }

  ## Hierarchy level 3: assign SNPs in introns/UTRs (excluding enhancer regions) to genes containing the introns/UTRs.
  if(!is.null(intron_utrs)){
    cat('Assign SNPs in UTRs (but not enhancer regions) to the UTR genes ... \n')
    if(intron.mode){
      cat('Assign SNPs in introns (but not in enhancer regions) to genes containing the introns ...\n')
    }
    intron_utr_overlap <- lapply(intron_utrs, function(x){plyranges::join_overlap_inner(x, finemap.gr)})
    snps.in <- unique(unlist(GenomicRanges::GRangesList(intron_utr_overlap))$snp)
    finemap.gr <- finemap.gr[!(finemap.gr$snp %in% snps.in),]
  }else{
    intron_utr_overlap <- NULL
  }

  ## Hierarchy level 4: assign the rest of SNPs (intergenic) to genes by distance weighting.
  if(length(finemap.gr) > 0){
    cat('Assign intergenic SNPs to genes by distance weighting ...\n')
    intergenic_snps_genes_by_distance <- list(intergenic = gene_by_distance(snps.gr = finemap.gr,
                                                                            promoters.gr = genomic.annots$promoters,
                                                                            c.dist = c.dist))
  }else{
    intergenic_snps_genes_by_distance <- NULL
  }

  snp.overlap.list <- c(exons_active_promoters_overlap,
                        enhancer_loops_overlap,
                        enhancer_snps_genes_by_distance,
                        intron_utr_overlap,
                        intergenic_snps_genes_by_distance)

  # Assign weights to SNP-gene pairs
  mat.list <- lapply(names(snp.overlap.list), function(x){
    snp.overlap.list[[x]] %>% as_tibble() %>%
      dplyr::mutate(category = x) %>%
      dplyr::distinct(gene_name, snp, category, .keep_all = T)
  })

  weights.mat <- Reduce(bind_rows, mat.list)
  weights.mat$weight <- ifelse( weights.mat$category %in% c(names(exons_active_promoters_overlap), names(enhancer_loops_overlap), names(intron_utr_overlap)),
                                1, weights.mat$weight )

  # Combine results from different categories
  weights.mat <- weights.mat %>%
    dplyr::arrange(category) %>%
    dplyr::group_by(snp, gene_name) %>%
    dplyr::mutate(category = paste(category, collapse = ',')) %>%
    dplyr::distinct(snp, gene_name, .keep_all = TRUE)

  # cat('SNP-gene pairs in each category: \n')
  # print(table(weights.mat$category))

  cat('Compute gene PIP ... \n')
  # For each SNP, distribute PIP of a SNP to all linked genes.
  # So weights of a SNP get normalized across genes (each SNP's weights of all genes should sum to 1)
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

  snp.gene.pip.mat <- suppressMessages(normalized.weights.mat %>% dplyr::left_join(., gene.pip.mat))
  snp.gene.pip.mat <- snp.gene.pip.mat[!is.na(snp.gene.pip.mat$gene_name),]

  snp.gene.pip.mat <- snp.gene.pip.mat %>%
    dplyr::select(all_of(cols.to.keep)) %>%
    as.data.frame()

  return(snp.gene.pip.mat)
}


#' @title Extract gene-level result from SNP-level gene mapping result
#'
#' @param snp.gene.pip.mat A data frame of SNP-level gene mapping result
#' @param gene.annots a GRanges object of gene annotations
#' @importFrom magrittr %>%
#' @return a data frame of gene-level view of gene mapping result
#' @export
extract_gene_level_result <- function(snp.gene.pip.mat, gene.annots) {
  cat('Extract gene level result ...\n')

  snp.gene.pip.mat <- snp.gene.pip.mat[!is.na(snp.gene.pip.mat$gene_name),]

  gene.locations <- as.data.frame(gene.annots)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]
  gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))

  m <- match(snp.gene.pip.mat$gene_name, gene.locations$gene_name)
  snp.gene.pip.mat$gene_chr <- gene.locations[m, 'seqnames']
  snp.gene.pip.mat$gene_pos <- gene.locations[m, 'tss']

  gene.pip.res <- snp.gene.pip.mat[, c('gene_chr', 'gene_pos', 'gene_name', 'gene_pip')] %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::rename(chr = gene_chr, pos = gene_pos) %>%
    as.data.frame()

  return(gene.pip.res)

}

#' @title Get credible gene sets from SNP level gene mapping table
#'
#' @param snp.gene.pip.mat A data frame of SNP level gene mapping table
#' @param by.locus Logical, if TRUE, get credible gene sets based on locus-level gene PIP,
#' If FALSE, get credible gene sets based on gene PIP.
#' @param gene.cs.percent.thresh percentage threshold for credible gene sets
#' @importFrom magrittr %>%
#' @return a data frame of credible gene set result. Columns:
#' gene_cs: credible gene sets,
#' gene_cs_locus_pip: credible gene sets and corresponding locus-level gene PIPs.
#' gene_cs_pip: credible gene sets and corresponding gene PIPs.
#' top_gene:  genes with highest gene PIP at each locus.
#' top_locus_gene_pip: locus-level gene PIP for the top gene.
#' Locus-level gene PIP only includes SNPs within a locus, so this value may be lower than the gene PIP.
#' top_gene_pip:  gene PIP of the top gene.
#' @export
gene_cs <- function(snp.gene.pip.mat,
                    by.locus = TRUE,
                    gene.cs.percent.thresh = 0.8){

  # Get locus level gene PIP
  # For each locus - gene pair, sum over the fractional PIPs for SNPs in the locus and linked to the gene
  snp.locus.gene.pip.mat <- snp.gene.pip.mat %>%
    dplyr::group_by(locus, gene_name) %>%
    dplyr::mutate(locus_gene_pip = sum(pip * frac_pip)) %>% dplyr::ungroup()

  # simplify to get locus, gene_name, locus_gene_pip, and gene_pip
  locus.gene.pip.df <- snp.locus.gene.pip.mat %>%
    dplyr::select(locus, gene_name, gene_pip, locus_gene_pip) %>%
    dplyr::distinct(locus, gene_name, .keep_all=TRUE)

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
      dplyr::arrange(desc(locus_gene_pip)) %>%
      dplyr::mutate(gene_pip_csum = cumsum(locus_gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.percent.thresh)[1])

    gene.cs.df <- gene.cumsum.df %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(gene_cs = paste0(gene_name, collapse=','),
                       gene_cs_locus_pip = paste(paste0(gene_name, '(',round(locus_gene_pip,3),')'), collapse=','),
                       top_gene = gene_name[1],
                       top_locus_gene_pip = locus_gene_pip[1],
                       top_gene_pip = gene_pip[1])
  }else{
    # for each locus, keep the genes with gene PIP cumsum > 0.8
    gene.cumsum.df <- locus.gene.pip.df %>%
      dplyr::group_by(locus) %>%
      dplyr::arrange(desc(gene_pip)) %>%
      dplyr::mutate(gene_pip_csum = cumsum(gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.percent.thresh)[1])

    gene.cs.df <- gene.cumsum.df %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(gene_cs = paste0(gene_name, collapse=','),
                       gene_cs_pip = paste(paste0(gene_name, '(',round(gene_pip,3),')'), collapse=','),
                       top_gene = gene_name[1],
                       top_gene_pip = gene_pip[1])
  }

  return(gene.cs.df)
}



#' @title Find all genes (promoters) within 1MB genes and assign weights based on distance
#'
#' @param snps.gr a GRanges object with the SNP locations.
#' @param promoters.gr a GRanges object with the gene promoter locations.
#' @param c.dist A scaling number used for computing weight based on SNP-gene distance.
#' Weight = exp(-dist/c). Default = 50000 (50kb).
#' @return a GRanges object with SNP to gene distance and weights calculated based on distance.
#' @export
#'
gene_by_distance <- function(snps.gr, promoters.gr, c.dist = 50000){

  if(is.null(promoters.gr$tss)){
    plus_strand <- which(as.character(strand(promoters.gr)) == '+')
    minus_strand <- which(as.character(strand(promoters.gr)) == '-')
    promoters.gr$tss <- NA
    promoters.gr$tss[plus_strand] <- end(promoters.gr)[plus_strand]
    promoters.gr$tss[minus_strand] <- start(promoters.gr)[minus_strand]
  }

  # get all promoters within 1MB
  res <- plyranges::join_overlap_inner(x = promoters.gr,
                                       y = snps.gr,
                                       maxgap = 1e6)

  res$distance <- abs(res$pos - res$tss)
  res$weight <- exp(-res$distance / as.numeric(c.dist))

  return(res)
}


#' @title Find the nearest genes for top SNPs in each locus.
#'
#' @param gwas.gr a GRanges object of the GWAS summary statistics
#' @param genes.gr a GRanges object of gene information
#' @param dist.to Find nearest genes by distance to gene body or TSS.
#' @param cols.to.keep columns to keep in the returned GRanges object
#' @importFrom magrittr %>%
#' @export
#' @return a data frame with SNP location and nearest gene.
#'
find_nearest_genes <- function(gwas.gr,
                               genes.gr,
                               dist.to = c('genebody', 'tss'),
                               cols.to.keep = c('snp','chr','pos', 'nearest_gene')){

  dist.to <- match.arg(dist.to)

  GenomeInfoDb::seqlevelsStyle(gwas.gr) <- 'UCSC'
  GenomeInfoDb::seqlevelsStyle(genes.gr) <- 'UCSC'

  if(!is.null(gwas.gr$zscore)){
    gwas.gr <- gwas.gr[order(abs(gwas.gr$zscore), decreasing = TRUE), ]
  }else{
    gwas.gr <- gwas.gr[order(gwas.gr$pval), ]
  }

  snps.gr <- gwas.gr[!duplicated(gwas.gr$locus), ]

  gene.locations <- as.data.frame(genes.gr)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]

  snps.gr$nearest_gene <- NA

  if(dist.to == 'tss'){
    gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))
    gene.locations$start <- gene.locations$end <- gene.locations$tss
    gene.locations.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations, keep.extra.columns = T)

    snp.nearest.gene.idx <- GenomicRanges::nearest(snps.gr, gene.locations.gr)
    snps.gr$nearest_gene <- gene.locations.gr$gene_name[snp.nearest.gene.idx]
  }else if(dist.to == 'genebody'){
    gene.locations.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations, keep.extra.columns = T)
    snp.nearest.gene.hits <- as.data.frame(nearest(snps.gr, gene.locations.gr, select = 'all'))
    colnames(snp.nearest.gene.hits) <- c('snp_idx', 'gene_idx')
    snp.nearest.gene.hits$snp <- snps.gr$snp[snp.nearest.gene.hits$snp_idx]
    snp.nearest.gene.hits$gene_name <- gene.locations.gr$gene_name[snp.nearest.gene.hits$gene_idx]

    snp.nearest.gene.df <- snp.nearest.gene.hits %>%
      dplyr::select(snp, gene_name) %>%
      dplyr::group_by(snp) %>%
      dplyr::summarise(nearestGene = paste(gene_name, collapse = ','))

    snp.nearest.gene.df$nearestGene <- as.character(snp.nearest.gene.df$nearestGene)

    snps.gr$nearest_gene <- snp.nearest.gene.df$nearestGene[match(snps.gr$snp, snp.nearest.gene.df$snp)]

  }

  res <- as.data.frame(snps.gr)[, cols.to.keep]
  return(res)
}


#' @title Gene view summary table
#'
#' @param genemapping_res data frame of gene mapping result
#' @param gene.pip.thresh Filter genes with gene PIP cutoff (default: 0.1)
#' @importFrom magrittr %>%
#' @export
#'
gene_view_summary <- function(genemapping_res, gene.pip.thresh = 0.1){
  gene.view.df <- genemapping_res %>%
    dplyr::mutate(fractional_PIP = pip * frac_pip) %>%
    dplyr::select(gene_name, gene_pip, fractional_PIP) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(gene_pip = round(gene_pip[1], 3),
                     n_snps_frac_pip_2percent = sum(fractional_PIP > 0.02)) %>%
    dplyr::filter(gene_pip > gene.pip.thresh)
  return(gene.view.df)
}

#' @title SNP view summary table
#'
#' @param genemapping_res data frame of gene mapping result
#' @param gene.annots data frame of gene annotations
#' @param finemap.gr GRange object of fine mapping result
#' @param fractional.PIP.thresh Filter SNPs with fractional PIP cutoff (default: 0.02)
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @export
#'
snp_view_summary <- function(genemapping_res, gene.annots, finemap.gr, fractional.PIP.thresh = 0.02){
  high.conf.snp.df <- genemapping_res %>% dplyr::filter(fractional_PIP > fractional.PIP.thresh)

  snp.gene <- high.conf.snp.df %>% dplyr::select(snp, pos, gene_name)

  gene.locs.df <- gene.annots %>% as_tibble()
  gene.locs.df$TSS <- ifelse(gene.locs.df$strand=='+', gene.locs.df$start, gene.locs.df$end)

  snp.gene.dist <- snp.gene %>%
    dplyr::left_join(., gene.locs.df, on = 'gene_name') %>%
    dplyr::mutate(dist = abs(TSS - pos)) %>%
    dplyr::select(snp, gene_name, dist)

  high.conf.snp.df <- dplyr::inner_join(high.conf.snp.df, snp.gene.dist, on = c('snp','gene_name'))

  # Add nearest gene (distance to gene body)
  GenomeInfoDb::seqlevelsStyle(finemap.gr) <- 'UCSC'
  snp.nearest.gene.gr <- find_nearest_genes(finemap.gr, gene.annots, dist.to = 'genebody')
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
#' @param genemapping_res data frame of gene mapping result
#' @param finemap.gr GRange object of fine mapping result
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @export
#'
LDblock_view_summary <- function(genemapping_res, finemap.gr){

  # Gene CS based on locus level gene PIP
  gene.cs.l <- gene_cs(genemapping_res, by.locus = TRUE, gene.cs.percent.thresh = 0.8)
  gene.cs.df <- gene.cs.l$gene.cs.df
  gene.cumsum.df <- gene.cs.l$gene.cumsum.df
  locus.gene.pip.df <- gene.cs.l$locus.gene.pip.df

  # Add nearest genes to LD blocks
  GenomeInfoDb::seqlevelsStyle(finemap.gr) <- 'UCSC'
  finemap.gr <- finemap.gr[order(abs(finemap.gr$zscore), decreasing = T), ]
  top.snps.gr <- finemap.gr[!duplicated(finemap.gr$locus), ]
  nearest_genebody_genes.gr <- find_nearest_genes(top.snps.gr, gene.annots, dist.to = 'genebody')

  locus_topsnp_nearest_genes.df <- nearest_genebody_genes.gr %>% as_tibble() %>% dplyr::select(locus, nearest_gene)

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

