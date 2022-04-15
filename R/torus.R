#' @title Prepare Torus input files
#' @description Prepares two files necessary for running Torus:
#' z-score file and annotation file.
#' @param sumstats cleaned summary statistics
#' @param annotation_bed_files annotation files in BED format.
#' The bed file must have three columns: chr, start, end.
#' Chromosomes should be numeric (no "chr")
#' and should be in hg19/b37 format.
#' You will get wrong results if you use hg38 or other
#' (The reference panel we used is in hg19).
#' @param torus_annot_file Annotation file name.
#' @param torus_zscore_file z-score file name.
#' @return a list containing paths to the z-score file and annotation file.
#' @export
#' @examples
#' torus.files <- prepare_torus_input_files(sumstats, annotation_bed_files)
prepare_torus_input_files <- function(sumstats, annotation_bed_files,
                                      torus_annot_file='torus_annotations.txt.gz',
                                      torus_zscore_file='torus_zscore.txt.gz'){

  stopifnot(all(file.exists(annotation_bed_files)))

  required.cols <- c('snp','chr', 'pos', 'locus', 'zscore')
  if(!all( required.cols %in% colnames(sumstats))){
    stop(sprintf('Column: %s cannot be found in the summary statistics!',
                 required.cols[which(!required.cols %in% colnames(sumstats))]))
  }

  cat('Annotating SNPs...\n')
  snps.annots <- annotate_snps_binary(sumstats, annotations = annotation_bed_files, keep.annot.only=T)

  cat('Generating Torus input files ...\n')
  if(missing(torus_annot_file)){
    torus_annot_file <- tempfile(pattern = "torus_annotations", fileext = "txt.gz")
  }
  readr::write_tsv(snps.annots, file=torus_annot_file, col_names = T)
  cat('Wrote Torus annotation files to', torus_annot_file, '\n')

  if(missing(torus_zscore_file)){
    torus_zscore_file <- tempfile(pattern = "torus_zscore", fileext = "txt.gz")
  }
  readr::write_tsv(sumstats[,c('snp','locus','zscore')], file=torus_zscore_file, col_names = T)
  cat('Wrote Torus z-score files to', torus_zscore_file, '\n')
  return(list(torus_annot_file=torus_annot_file, torus_zscore_file=torus_zscore_file))
}

#' @title Run enrichment analysis and compute SNP-level priors using Torus
#' @description Perform enrichment analysis using Torus and then
#' compute SNP-level priors using the enrichment estimates.
#' @param torus_annot_file SNP annotation file prepared by
#' the \code{prepare_torus_input_files} function.
#' The SNP annotation file contains SNP-level genomic annotations used by
#' Torus analysis. The annotation file uses a header to specify the number
#' and the nature (categorical or continuous) of the annotations.
#' The first column with the header "SNP" represents the SNP name.
#' The following columns represent specific annotations.
#' For categorical/discrete annotations, the header should have a suffix "_d";
#' whereas for continuous annotations, the header should ends with "_c".
#' @param torus_zscore_file Summary statistics from single SNP association
#' analysis, prepared by
#' the \code{prepare_torus_input_files} function.
#' Should be compressed in gzip format.
#' @param option Torus options:
#' \dQuote{est}, obtain estimates of enrichment parameters and their confidence intervals;
#' \dQuote{est-prior}, perform enrichment analysis and
#' compute SNP-level priors using
#' the estimated enrichment estimates for each locus;
#' or \dQuote{fdr}, perform Bayesian FDR control, and output the result.
#' @param torus_path Path to \code{torus} executable.
#' @importFrom tibble as_tibble
#' @return a list of enrichment results and SNP-level prior probabilities.
#' Enrichment result contains the point estimate (MLE) of the log odds ratio,
#' as well as 95% confidence interval for the corresponding point estimate.
#' @export
#' @examples
#' # Get enrichment estimates and confidence intervals
#' torus.result <- run_torus("torus_annotations.txt.gz",
#'                           "torus_zscore.txt.gz",
#'                           option = "est")
#'
#' # Get enrichment estimates and compute SNP-level priors
#' torus.result <- run_torus("torus_annotations.txt.gz",
#'                           "torus_zscore.txt.gz",
#'                           option = "est-prior")
#' # Bayesian FDR control
#' torus.result <- run_torus("torus_annotations.txt.gz",
#'                           "torus_zscore.txt.gz",
#'                           option = "fdr")
#'
run_torus <- function(torus_annot_file,
                      torus_zscore_file,
                      option=c('est', 'est-prior', 'fdr'),
                      torus_path='torus'){

  if(!file.exists(torus_annot_file)){
    stop('Cannot find Torus annotation file.')
  }
  if(!file.exists(torus_zscore_file)){
    stop('Cannot find Torus z-score file.')
  }

  option <- match.arg(option)
  cat('Run Torus...\n')

  torus.result <- list()
  if(option == 'est'){
    torus_args <- c('-d', torus_zscore_file,
                    '-annot', torus_annot_file,
                    '--load_zval',
                    '-est')
    cat('Estimating enrichments...\n')
    res <- processx::run(command = torus_path, args = torus_args, echo_cmd = TRUE, echo = TRUE)
    enrich <- as_tibble(read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F))
    colnames(enrich) <- c("term", "estimate", "low", "high")
    return(enrich)

  }else if(option == 'est-prior'){
    torus_args <- c('-d', torus_zscore_file,
                    '-annot', torus_annot_file,
                    '--load_zval',
                    '-est',
                    '-dump_prior', 'prior')
    cat('Estimating enrichments and computing SNP-level priors...\n')
    res <- processx::run(command = torus_path, args = torus_args, echo_cmd = TRUE, echo = TRUE)
    enrich <- as_tibble(read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F))
    colnames(enrich) <- c("term", "estimate", "low", "high")

    files <- list.files(path = 'prior/', pattern = '*.prior', full.names = T)
    files.str <- paste0(files, collapse = " ")
    system(paste('cat', files.str, '> prior/allchunks.txt'))
    snp_pip <- suppressMessages(vroom::vroom('prior/allchunks.txt', col_names = F, delim = "  "))
    colnames(snp_pip) <- c("snp","torus_pip")
    system('rm -rf prior/')
    return(list(enrich = enrich, snp_pip = snp_pip))

  }else if(option == 'fdr'){
    torus_args <- c('-d', torus_zscore_file,
                    '-annot', torus_annot_file,
                    '--load_zval',
                    '-qtl')
    cat('Performing Bayesian FDR control...\n')
    res <- processx::run(command = torus_path, args = torus_args, echo_cmd = TRUE, echo = TRUE)
    torus_fdr <- as_tibble(read.table(file = textConnection(res$stdout),header=F,stringsAsFactors = F))
    colnames(torus_fdr) <- c("rej","region_id","fdr","decision")
    return(torus_fdr)
  }
}

