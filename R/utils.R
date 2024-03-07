
#' @title Process PC-HiC data and save as a GRanges object
#'
#' @param pcHiC A data frame of PC-HiC data,
#' with columns of enhancer"Promoter" and
#' "Interacting_fragment". Interacting_fragment should contains
#' chr, start and end positions of the fragments interacting with promoters
#' e.g. "chr.start.end" or "chr:start-end".
#' @param gene.annots If provided, restrict to genes in gene.annots
#' @param score.thresh Numeric. Threshold of interaction scores.
#' (default = 0).
#' @param flank  Integer. Extend bases on both sides of the regulatory elements
#' (default = 0).
#' @return A GRanges object with processed PC-HiC links, with genomic coordinates
#' of the interacting regions and gene names (promoters).
#' @export
process_pcHiC <- function(pcHiC, gene.annots, score.thresh = 0, flank = 0){

  pcHiC <- pcHiC %>% dplyr::select(Promoter, Interacting_fragment)
  # separate genes connecting to the same fragment
  pcHiC <- pcHiC %>%
    tidyr::separate_rows(Promoter) %>%
    dplyr::rename(gene_name = Promoter)

  pcHiC <- pcHiC %>%
    tidyr::separate(Interacting_fragment, c('chr', 'start', 'end')) %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))

  if(!missing(gene.annots)){
    pcHiC <- pcHiC %>% dplyr::filter(gene_name %in% gene.annots$gene_name)
  }

  if(score.thresh > 0){
    pcHiC <- pcHiC %>% dplyr::filter(score >= score.thresh)
  }

  if(flank > 0){
    pcHiC <- pcHiC %>% dplyr::mutate(start = start - as.integer(flank),
                                     end = end + as.integer(flank))
  }

  columns <- c('chr', 'start', 'end', 'promoter_start', 'promoter_end', 'gene_name', 'score')
  pcHiC <- pcHiC %>% dplyr::select(columns)

  pcHiC.gr <- GenomicRanges::makeGRangesFromDataFrame(pcHiC, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(pcHiC.gr) <- 'UCSC'

  return(pcHiC.gr)
}


#' @title Process ABC scores and save as a GRanges object
#'
#' @param ABC A data frame of ABC scores from Nasser et al. Nature 2021 paper
#' @param gene.annots If provided, restrict to genes in gene.annots
#' @param full.element Logical; if TRUE, use full length of ABC elements
#' extracted from the "name" column. Otherwise, use the original (narrow)
#' regions provided in the ABC scores data.
#' @param score.thresh Numeric. Threshold of ABC scores.
#' (default = 0.015, as in Nasser et al. Nature 2021 paper).
#' @param flank  Integer. Extend bases on both sides of the ABC elements (default = 0).
#' @return a GRanges object with processed ABC scores, with genomic coordinates
#' of the interacting regions and gene names (promoters).
#' @export
process_ABC <- function(ABC, gene.annots, full.element = FALSE, score.thresh = 0.015, flank = 0){

  if(full.element){
    ABC <- ABC %>%
      tidyr::separate(name, c(NA, 'element_region'), sep = '\\|', remove = FALSE) %>%
      tidyr::separate(element_region, c(NA, 'element_location'), sep = '\\:') %>%
      tidyr::separate(element_location, c('element_start', 'element_end'), sep = '\\-') %>%
      dplyr::mutate(start = as.numeric(element_start), end = as.numeric(element_end),
                    promoter_start = TargetGeneTSS,
                    promoter_end = TargetGeneTSS)
  }

  ABC <- ABC %>% dplyr::rename(gene_name = TargetGene, score = ABC.Score)

  if(!missing(gene.annots)){
    ABC <- ABC %>% dplyr::filter(gene_name %in% gene.annots$gene_name)
  }

  if(score.thresh > 0){
    ABC <- ABC %>% dplyr::filter(score >= score.thresh)
  }

  if(flank > 0){
    ABC <- ABC %>% dplyr::mutate(start = start - as.integer(flank),
                                 end = end + as.integer(flank))
  }

  columns <- c('chr', 'start', 'end', 'promoter_start', 'promoter_end', 'gene_name', 'score')
  ABC <- ABC %>% dplyr::select(columns)

  ABC.gr <- GenomicRanges::makeGRangesFromDataFrame(ABC, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(ABC.gr) <- 'UCSC'

  return(ABC.gr)
}

#' @title Process chromatin loop data and save as a GRanges object
#'
#' @param loops A data frame of chromatin loops, with columns:
#' "chr", "start", "end", "gene_name", and "score" (optional).
#' @param gene.annots If provided, restrict to genes in gene.annots
#' @param score.thresh Numeric. Threshold of interaction scores.
#' (default = 0).
#' @param flank  Integer. Extend bases on both sides of the regulatory elements
#' (default = 0).
#' @return A GRanges object with processed chromatin loops,
#' with genomic coordinates of the regulatory elements and gene names.
#' @export
process_loop_data <- function(loops, gene.annots, score.thresh = 0, flank = 0){

  loops <- as.data.frame(loops)

  if(!missing(gene.annots)){
    loops <- loops %>% dplyr::filter(gene_name %in% gene.annots$gene_name)
  }

  if(score.thresh > 0){
    loops <- loops %>% dplyr::filter(score >= score.thresh)
  }

  if(flank > 0){
    loops <- loops %>% dplyr::mutate(start = start - as.integer(flank),
                                     end = end + as.integer(flank))
  }

  loops.gr <- GenomicRanges::makeGRangesFromDataFrame(loops, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(loops.gr) <- 'UCSC'

  return(loops.gr)
}


#' @title Get nearby interactions for enhancer regions near promoters
#'
#' @param enhancer_regions A GRanges object of enhancer regions
#' @param promoters A GRanges object of promoters
#' @param max.dist Max distance between enhancer regions and promoters (default: 20kb)
#' @return A GRanges object of nearby interactions
#' @export
nearby_interactions <- function(enhancer_regions, promoters, max.dist = 20000){

  promoters$promoter_chr <- seqnames(promoters)
  promoters$promoter_start <- start(promoters)
  promoters$promoter_end <- end(promoters)

  nearby.interactions <- plyranges::join_overlap_inner(enhancer_regions,
                                                       promoters,
                                                       maxgap = max.dist)

  columns <- c('promoter_chr', 'promoter_start', 'promoter_end', 'gene_name')
  nearby.interactions <- nearby.interactions[, columns]

  return(nearby.interactions)
}

# Load ENCODE narrow peak data and convert to a GRanges object
process_narrowpeaks <- function(peak.file){
  if( !file.exists(peak.file) ){stop('narrowpeak file is not availble!')}
  peaks <- data.table::fread(peak.file)
  colnames(peaks) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevelsStyle(peaks.gr) <- 'UCSC'
  return(peaks.gr)
}


# Annotations for causal SNPs
annotator_merged <- function(sumstats, annotations){

  snpRanges <- make_ranges(sumstats$chr, sumstats$pos, sumstats$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=sumstats$snp)
  sumstats['annots'] <- ''

  for(f in annotations){

    curr <- rtracklayer::import(f, format='bed')
    snpRangesIn <- IRanges::subsetByOverlaps(snpRanges, annot.gr)
    snpsIn <- unique(snpRangesIn$snp)

    if(length(snpsIn)>0){
      curr <- sumstats %>% dplyr::pull(annots)
      curr <- curr[sumstats$snp %in% snpsIn]
      delims <- rep(';', length(curr))
      delims[which(curr == '')] <- ''
      sumstats[sumstats$snp %in% snpsIn,'annots'] <-
        paste0(curr,delims,gsub(pattern = '.bed',replacement = '', x = basename(f)))
    }
  }
  return(sumstats)
}

# helper function to make a GRanges object
make_ranges <- function(seqname, start, end){
  return(GenomicRanges::GRanges(seqnames = seqname,
                                ranges = IRanges::IRanges(start = start, end = end)))
}

# get TSS positions from promoters or gene annotations.
get_tss <- function(gr, type = c('promoter', 'gene')){
  type <- match.arg(type)

  if(type == 'gene'){
    tss <- GenomicRanges::start(GenomicRanges::resize(gr, width = 1))
  }else if(type == 'promoter'){
    tss <- integer(length(gr))
    plus_strand <- which(as.character(GenomicRanges::strand(gr)) == '+')
    minus_strand <- which(as.character(GenomicRanges::strand(gr)) == '-')
    tss[plus_strand] <- GenomicRanges::end(gr)[plus_strand]
    tss[minus_strand] <- GenomicRanges::start(gr)[minus_strand]
  }

  return(tss)
}

# get LD (r^2) between each SNP and the top SNP in sumstats from bigSNP
get_LD_bigSNP <- function(sumstats, bigSNP, topSNP = NULL){

  # only include SNPs in bigSNP markers
  sumstats <- sumstats[sumstats$snp %in% bigSNP$map$marker.ID, ]
  sumstats$bigSNP_idx <- match(sumstats$snp, bigSNP$map$marker.ID)

  if(missing(topSNP)){
    if( max(sumstats$pval) <= 1 ){
      sumstats$pval <- -log10(sumstats$pval)
    }
    top_snp_idx <- sumstats$bigSNP_idx[which.max(sumstats$pval)]
  }else{
    top_snp_idx <- sumstats$bigSNP_idx[sumstats$snp == topSNP]
  }

  top_snp_genotype <- bigSNP$genotypes[,top_snp_idx]
  genotype.mat <- bigSNP$genotypes[,sumstats$bigSNP_idx]

  r2.vals <- as.vector(cor(top_snp_genotype, genotype.mat))^2
  sumstats$r2 <- round(r2.vals, 4)

  return(sumstats)
}


# Add LD information from bigSNP
add_LD_bigSNP <- function(sumstats, bigSNP,
                          r2.breaks = c(0, 0.1, 0.25, 0.75, 0.9, 1),
                          r2.labels = c('0-0.1','0.1-0.25','0.25-0.75','0.75-0.9','0.9-1')) {

  # only include SNPs in bigSNP markers
  sumstats <- sumstats[sumstats$snp %in% bigSNP$map$marker.ID, ]
  sumstats$bigSNP_idx <- match(sumstats$snp, bigSNP$map$marker.ID)

  if( max(sumstats$pval) <= 1 ){
    sumstats$pval <- -log10(sumstats$pval)
  }

  locus_list <- unique(sumstats$locus)

  sumstats.r2.df <- data.frame()
  for(locus in locus_list){
    curr_sumstats <- sumstats[sumstats$locus == locus, ]
    top_snp_idx <- curr_sumstats$bigSNP_idx[which.max(curr_sumstats$pval)]
    top_snp_genotype <- bigSNP$genotypes[,top_snp_idx]
    genotype.mat <- bigSNP$genotypes[,curr_sumstats$bigSNP_idx]

    r2.vals <- as.vector(cor(top_snp_genotype, genotype.mat))^2
    r2.brackets <- cut(r2.vals, breaks = r2.breaks, labels = r2.labels)
    curr_sumstats$r2 <- r2.brackets
    sumstats.r2.df <- rbind(sumstats.r2.df, curr_sumstats)
  }

  return(sumstats.r2.df)
}

# get gene region
get_gene_region <- function(gene.mapping.res,
                            genes.of.interest,
                            ext = 10000,
                            select.region = c('all', 'locus')){

  select.region <- match.arg(select.region)

  high.conf.snp.df <- gene.mapping.res %>% dplyr::filter(pip > 0.2) %>%
    dplyr::group_by(snp) %>% dplyr::arrange(-gene_pip) %>% dplyr::slice(1)
  gene.gr <- gene.annots[match(high.conf.snp.df$gene_name, gene.annots$gene_name),]
  gene.gr$tss <- start(resize(gene.gr, width = 1))
  gene.gr <- gene.gr[,c('gene_name','tss')]
  high.conf.snp.df$tss <- gene.gr$tss

  gene.snp.tss <- high.conf.snp.df %>%
    dplyr::filter(gene_name %in% genes.of.interest) %>%
    dplyr::group_by(locus) %>%
    dplyr::arrange(-pip) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(distToTSS = pos-tss) %>%
    dplyr::select(gene_name, locus, chr, pos, tss, distToTSS)

  if(select.region == 'all'){
    chr <- gene.snp.tss$chr[1]
    locus <- paste(gene.snp.tss$locus, collapse = ',')
    region_start <- min(c(gene.snp.tss$pos, gene.snp.tss$tss)) - ext
    region_end <- max(c(gene.snp.tss$pos, gene.snp.tss$tss)) + ext
    region <- GRanges(seqnames = chr,
                      IRanges(start = region_start, end = region_end),
                      locus = locus)
    seqlevelsStyle(region) <- 'UCSC'
  }else if(select.region == 'locus'){
    region <- lapply(gene.snp.tss$locus, function(l){
      locus <- l
      locus.snp.tss <- gene.snp.tss[gene.snp.tss$locus == locus, ]
      chr <- locus.snp.tss$chr
      if(distToTSS < 0){ # snp is upstream
        region_start <- locus.snp.tss$pos - ext
        region_end <- locus.snp.tss$tss + ext
      } else{ # snp is downstream
        region_start <- locus.snp.tss$tss - ext
        region_end <- locus.snp.tss$pos + ext
      }
      gene_locus_region.gr <- GRanges(seqnames = chr,
                                      IRanges(start = region_start, end = region_end),
                                      locus = locus)
      seqlevelsStyle(gene_locus_region.gr) <- 'UCSC'
      gene_locus_region.gr
    })
  }
  return(region)
}

#' read LD matrix data by file format
read_LD <- function(file, format = c("rds", "rdata", "csv", "txt", "tsv")) {
  format <- match.arg(format)

  # if format is missing, try to guess format by file extension
  if (missing(format)) {
    file_ext_lower <- tolower(tools::file_ext(file))

    if (file_ext_lower == "rds"){
      format <- "rds"
    } else if (file_ext_lower %in% c("rdata", "rd", "rda", "rdat")){
      format <- "rd"
    } else if (file_ext_upper %in% c("csv", "csv.gz")) {
      format <- "csv"
    } else if (file_ext_lower %in% c("txt", "txt.gz")){
      format <- "txt"
    } else if (file_ext_lower %in% c("tsv", "tsv.gz")){
      format <- "tsv"
    } else {
      stop("Unknown file format!")
    }
  }

  if (format == "rds"){
    R <- readRDS(file)
  } else if (format == "rd"){
    R <- get(load(file))
  } else if (format == "csv"){
    R <- as.matrix(read.csv(file, sep=",", row.names=1))
  } else if (format %in% c("txt", "tsv")){
    R <- as.matrix(data.table::fread(file))
  } else {
    stop("Unknown file format!")
  }

  return(R)
}


#' Get region info with filenames of LD matrices and variant information
#'
#' @param LD_Blocks A data frame of LD blocks
#' @param LDREF.dir Directory of UKBB LD reference files
#' @param prefix prefix name of the UKBB LD reference files
#' @param LD_matrix_ext File extension of LD matrix files
#' @param snp_info_ext File extension of SNP information files
#' @return A data frame with information of the variants in the LD matrix.
#' @export
get_UKBB_region_info <- function(LD_Blocks,
                                 LDREF.dir,
                                 prefix = "ukb_b37_0.1",
                                 LD_matrix_ext = "RDS",
                                 snp_info_ext = "Rvar") {

  LD.file <- sprintf("%s_chr%d.R_snp.%d_%d", prefix, LD_Blocks$chr, LD_Blocks$start, LD_Blocks$end)
  LD_matrix_file <- file.path(LDREF.dir, paste0(LD.file, ".", LD_matrix_ext))
  snp_info_file <- file.path(LDREF.dir, paste0(LD.file, ".", snp_info_ext))
  region_info <- data.frame(LD_Blocks,
                            LD_matrix = LD_matrix_file,
                            snp_info = snp_info_file)

  return(region_info)
}

#' Load UK Biobank LD reference matrix and variant information
#'
#' @param LD_Blocks A data frame of LD blocks
#' @param locus locus ID
#' @param LDREF.dir Directory of UKBB LD reference files
#' @param prefix prefix name of the UKBB LD reference files
#' @param LD_matrix_ext File extension of LD matrix files
#' @param snp_info_ext File extension of SNP information files
#'
#' @return A list, containing LD (correlation) matrix R and
#' a data frame with information of the variants in the LD matrix.
#' @export
load_UKBB_LDREF <- function(LD_Blocks,
                            locus,
                            LDREF.dir,
                            prefix = "ukb_b37_0.1",
                            LD_matrix_ext = "RDS",
                            snp_info_ext = "Rvar"){
  if(!locus %in% LD_Blocks$locus){
    stop("locus is not in LD_blocks!")
  }
  LD_Block <- LD_Blocks[LD_Blocks$locus == locus, ]
  LD.file <- sprintf("%s_chr%d.R_snp.%d_%d", prefix, LD_Block$chr, LD_Block$start, LD_Block$end)
  R <- readRDS(file.path(LDREF.dir, paste0(LD.file, ".", LD_matrix_ext)))
  snp_info <- data.table::fread(file.path(LDREF.dir, paste0(LD.file, ".", snp_info_ext)))
  return(list('R' = R, 'snp_info' = snp_info))
}

