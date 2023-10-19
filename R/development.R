
# Subtracts all exon coordinates from all gene coordinates to get all intron coordinates
subtract_exons <- function( genes.gr, exons.gr ) {
  # adapted from bosberg on Biostars. (https://www.biostars.org/p/489350/)
  # Subtract GRange object gr2 from gr1, but unlike setdiff, preserve individual ranges in gr1
  genes.df <- data.frame( seqnames=seqnames(genes.gr),
                          start=start(genes.gr)-1, end=end(genes.gr), # convert to 0-based start for bedtools
                          strand=strand(genes.gr),
                          gene_name = genes.gr$gene_name )
  exons.df <- data.frame( seqnames=seqnames(exons.gr),
                          start=start(exons.gr)-1,
                          end=end(exons.gr),
                          strand=strand(exons.gr) )
  result <- bedtoolsr::bt.subtract(genes.df, exons.df)

  if ( length(result)==0 ){
    # subtraction has left nothing remaining. Return empty GRanges obj.
    return( GRanges() )
  } else {
    colnames(result) <- colnames(genes.df)
    result$start <- result$start+1 # reset to 1-based notation consistent with GRanges
    return( GRanges(result) )
  }
}
