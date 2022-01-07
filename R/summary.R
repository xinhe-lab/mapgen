
#' @title Gene view summary table
#'
#' @param genemapping_res data frame of gene mapping result
#' @param gene.pip.thresh Filter genes with gene PIP cutoff (default: 0.1)
#' @import tidyverse
#'
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
#' @import tidyverse
#'
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
#' @import tidyverse
#'
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
