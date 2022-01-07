#' @title PrepareTorusFiles
#' @description Prepares two files necessary for running torus: z-score file and annotations file
#' @param cleaned_sumstats cleaned summary statistics from RunCleaner or other
#' @param bed_annotation a directory containing annotations in bed format. The bed file must have three columns: chr, start, end. Chromosomes should be numeric (no "chr")
#' and should be in hg19/b37 format. You will get wrong results if you use hg38 or other (reference panel is hg19).
#' @return location to two files: torus zscore file and annotation file
#' @export
PrepareTorusFiles <- function(cleaned_sumstats, bed_annotations){
  
  stopifnot(dir.exists(bed_annotations))

  annotations <- list.files(path = bed_annotations, pattern = '*.bed', full.names = T)
  
  if(length(annotations) == 0){
    stop('No bed files/annotations were found in this directory. Are you sure this directory contains bed files?')
  }
  
  print('Annotating Snps..')
  cleaned.gwas.annots <- annotator(cleaned_sumstats, annotations = annotations)
  
  print('Writing files to temporary location..')
  
  readr::write_tsv(x = cleaned.gwas.annots[,-c(1:6,8:12)], path = 'torus_annotations.txt.gz', col_names = T)
  readr::write_tsv(x = cleaned_sumstats[,c('snp','locus','zscore')], path = 'torus_zscore.txt.gz', col_names = T)
  
  print('Done.')
  
  return(list(torus_annot_file='torus_annotations.txt.gz', torus_zscore_file='torus_zscore.txt.gz'))
}

#' @title RunTorus
#' @description executes Torus
#' @return list of enrichment results and PIP of each SNP (both tibbles)
#' @export
RunTorus <- function(torus_annot_file, torus_zscore_file, TORUS=system.file('torus', package='finemappeR')){
  
  if(!file.exists(torus_annot_file)){
    stop('Cannot find annotation files. Did you run PrepareTorusfiles?')
  }
  if(!file.exists(torus_zscore_file)){
    stop('Cannot find zscore files. Did you run PrepareTorusfiles?')
  }
  
  args <- c('-d',
            torus_zscore_file, 
            '-annot', 
            torus_annot_file,
            '--load_zval',
            '-dump_prior',
            'prior')
  
  res <- processx::run(command = TORUS, args = args, echo_cmd = TRUE, echo = TRUE)
  print('Saving enrichments to dataframe..')
  enrich <- as_tibble(read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F))
  colnames(enrich) <- c("term", "estimate", "low", "high")
  
  print('Extracting prior probabilities from Torus..')
  files <- list.files(path = 'prior/', pattern = '*.prior', full.names = T)
  files.str <- paste0(files, collapse = " ")
  
  system(paste0('cat ', files.str, ' > prior/allchunks.txt'))
  snp_pip <- suppressMessages(vroom::vroom('prior/allchunks.txt', col_names = F, delim = "  "))
  colnames(snp_pip) <- c("snp","torus_pip")
  
  system('rm -rf prior/')
  
  return(list(enrich=enrich, snp_pip=snp_pip))
  
}

#' @title RunTorusFDR
#' @description Runs Torus in FDR mode to estimate the FDR of each chunk containing a causal variant
#' @return tibble containing fdr of each LD block/chunk
#' @export
RunTorusFDR <- function(torus_annot_file, torus_zscore_file, TORUS=system.file('torus', package='finemappeR')){
  
  if(!dir.exists('.temp')){
    stop('Cannot find annotation files. Did you run PrepareTorusfiles?')
  }
  
  args <- c('-d',
            torus_zscore_file, 
            '-annot', 
            torus_annot_file,
            '--load_zval',
            '-qtl')
  
  res <- processx::run(command = TORUS, args = args, echo_cmd = TRUE, echo = TRUE)
  torus_fdr <- as_tibble(read.table(file = textConnection(res$stdout),header=F,stringsAsFactors = F))
  colnames(torus_fdr) <- c("rej","region_id","fdr","decision")

  return(torus_fdr)
}

