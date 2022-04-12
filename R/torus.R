#' @title Prepare Torus input files
#' @description Prepares two files necessary for running Torus:
#' z-score file and annotation file.
#' @param cleaned_sumstats cleaned summary statistics from RunCleaner or other
#' @param bed_annotation a directory containing annotations in bed format.
#' The bed file must have three columns: chr, start, end.
#' Chromosomes should be numeric (no "chr")
#' and should be in hg19/b37 format.
#' You will get wrong results if you use hg38 or other
#' (The reference panel we used is in hg19).
#' @param torus_annot_file Annotation file name.
#' @param torus_zscore_file z-score file name.
#' @return a list containing paths to the z-score file and annotation file.
#' @export
PrepareTorusFiles <- function(cleaned_sumstats, bed_annotations,
                              torus_annot_file='torus_annotations.txt.gz',
                              torus_zscore_file='torus_zscore.txt.gz'){

  stopifnot(dir.exists(bed_annotations))

  annotations <- list.files(path = bed_annotations, pattern = '*.bed', full.names = TRUE)

  if(length(annotations) == 0){
    stop('No bed files/annotations were found in this directory. Are you sure this directory contains bed files?')
  }

  cat('Annotating SNPs...\n')
  cleaned.gwas.annots <- annotator(cleaned_sumstats, annotations = annotations)

  cat('Writing Torus input files...\n')
  readr::write_tsv(cleaned.gwas.annots[,-c(1:6,8:12)], file=torus_annot_file, col_names = T)
  readr::write_tsv(cleaned_sumstats[,c('snp','locus','zscore')], file=torus_zscore_file, col_names = T)

  return(list(torus_annot_file=torus_annot_file, torus_zscore_file=torus_zscore_file))
}

#' @title Run enrichment analysis and compute SNP-level priors using Torus
#' @description Perform enrichment analysis using Torus and then
#' compute SNP-level priors using the enrichment estimates.
#' @param torus_annot_file SNP annotation file.
#' The SNP annotation file contains SNP-level genomic annotations used by
#' TORUS analysis. The annotation file uses a header to specify the number
#' and the nature (categorical or continuous) of the annotations.
#' The first column with the header "SNP" represents the SNP name.
#' The following columns represent specific annotations.
#' For categorical/discrete annotations, the header should have a suffix "_d";
#' whereas for continuous annotations, the header should ends with "_c".
#' @param torus_zscore_file Summary statistics from single SNP association
#' analysis. Should be compressed in gzip format.
#' @param TORUS Path to Torus executable.
#' @param dump_prior Logical. If TURE, compute SNP-level priors using
#' the estimated enrichment estimates for each locus.
#' @importFrom tibble as_tibble
#' @return a list of enrichment results and SNP-level prior probabilities.
#' Enrichment result contains the point estimate (MLE) of the log odds ratio,
#' as well as 95% confidence interval for the corresponding point estimate.
#' @export
RunTorus <- function(torus_annot_file,
                     torus_zscore_file,
                     TORUS='torus',
                     dump_prior=TRUE){

  if(!file.exists(torus_annot_file)){
    stop('Cannot find annotation files. Did you run PrepareTorusfiles?')
  }
  if(!file.exists(torus_zscore_file)){
    stop('Cannot find zscore files. Did you run PrepareTorusfiles?')
  }

  args <- c('-d', torus_zscore_file,
            '-annot', torus_annot_file,
            '--load_zval',
            '-est')

  if(dump_prior){
    args <- c(args, '-dump_prior', 'prior')
  }

  cat('Estimating enrichments...\n')
  res <- processx::run(command = TORUS, args = args, echo_cmd = TRUE, echo = TRUE)
  enrich <- as_tibble(read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F))
  colnames(enrich) <- c("term", "estimate", "low", "high")

  if(dump_prior){
    cat('Extracting prior probabilities from Torus...\n')
    files <- list.files(path = 'prior/', pattern = '*.prior', full.names = T)
    files.str <- paste0(files, collapse = " ")

    system(paste('cat', files.str, '> prior/allchunks.txt'))
    snp_pip <- suppressMessages(vroom::vroom('prior/allchunks.txt', col_names = F, delim = "  "))
    colnames(snp_pip) <- c("snp","torus_pip")
    system('rm -rf prior/')
  }else{
    snp_pip <- NULL
  }

  return(list(enrich=enrich, snp_pip=NULL))

}

#' @title Perform Bayesian FDR control using Torus
#' @description Runs Torus in FDR mode to estimate the FDR of
#' each chunk containing a causal variant.
#' @param torus_annot_file SNP annotation file.
#' The SNP annotation file contains SNP-level genomic annotations used by
#' TORUS analysis. The annotation file uses a header to specify the number
#' and the nature (categorical or continuous) of the annotations.
#' The first column with the header "SNP" represents the SNP name.
#' The following columns represent specific annotations.
#' For categorical/discrete annotations, the header should have a suffix "_d";
#' whereas for continuous annotations, the header should ends with "_c".
#' @param torus_zscore_file Summary statistics from single SNP association
#' analysis. Should be compressed in gzip format.
#' @param TORUS Path to Torus executable.
#' @importFrom tibble as_tibble
#' @return a tibble containing FDR of each LD block/chunk
#' @export
RunTorusFDR <- function(torus_annot_file, torus_zscore_file, TORUS='torus'){

  if(!dir.exists('.temp')){
    stop('Cannot find annotation files. Did you run PrepareTorusfiles?')
  }

  args <- c('-d', torus_zscore_file,
            '-annot', torus_annot_file,
            '--load_zval',
            '-qtl')

  res <- processx::run(command = TORUS, args = args, echo_cmd = TRUE, echo = TRUE)
  torus_fdr <- as_tibble(read.table(file = textConnection(res$stdout),header=F,stringsAsFactors = F))
  colnames(torus_fdr) <- c("rej","region_id","fdr","decision")

  return(torus_fdr)
}

