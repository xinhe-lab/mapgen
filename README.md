
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Mapgen

<!-- badges: start -->
<!-- badges: end -->

Mapgen is a multi-function software that performs the following tasks:

1.  enrichment analysis of functional annotations for a trait of
    interest.
2.  functionally-informed genetic fine-mapping.
3.  gene mapping based on fine-mapping result and genomic annotations.

## Installation

You can install the development version of `mapgen` from
[GitHub](https://github.com/xinhe-lab/mapgen) with:

``` r
install.packages("remotes")
remotes::install_github("xinhe-lab/mapgen")
```

After installing, check that it loads properly:

``` r
library(mapgen)
```

## Main steps

### 1. Data preparation:

Prepare GWAS summary statistics, functional annotations, as well as LD
reference panel as input data.

### 2. Enrichment analysis:

Assess the enrichment of genetic signals of a trait of interest in
functional annotations.

\*Please install [TORUS](https://github.com/xqwen/torus) software
package, if you need to run enrichment analysis.

### 3. Fine-mapping:

Perform Bayesian statistical fine-mapping using SuSiE on
trait-associated loci, using a informative prior that favors variants
located in enriched annotations.

\*Please install [susieR](https://github.com/stephenslab/susieR)
package, if you need to run fine-mapping with GWAS summary statistics
using SuSiE.

### 4. Gene mapping:

Infer causal genes at each locus based on fine-mapping result and
genomic annotations, including gene annotations, chromatin loops, etc.

## Example tutorials using data from our AFib study

1.  Data preparation: obtain AFib GWAS data and cell-type-resolved open
    chromatin regions (OCRs) from scATAC-seq.
2.  [Enrichment
    analysis](https://xinhe-lab.github.io/mapgen/articles/enrichment_finemapping_tutorial.html):
    estimate the enrichment of AFib signals in OCRs across cell types.
3.  [Fine-mapping](https://xinhe-lab.github.io/mapgen/articles/enrichment_finemapping_tutorial.html):
    Perform fine-mapping on AFib-associated loci, using a informative
    prior that favors variants located in OCRs of enriched cell types.
4.  [Partition fine-mapping PIPs by annotation
    categories](https://xinhe-lab.github.io/mapgen/articles/partition_pip_tutorial.html):
    Assign the likely cell type(s) through which the causal variants act
    in most loci using fine-mapped SNPs and its associated cell type
    information.
5.  [Gene
    mapping](https://xinhe-lab.github.io/mapgen/articles/gene_mapping_tutorial.html):
    infer causal genes (gene PIPs) at each locus based on AFib
    fine-mapping result and genomic annotations, including gene
    annotations, chromatin loops (PC-HiC links, ABC scores), etc.

<img src="man/figures/workflow.overview.png" width="75%" />

## Reference

Alan Selewa\*, Kaixuan Luo\*, Michael Wasney, Linsin Smith, Xiaotong
Sun, Chenwei Tang, Heather Eckart, Ivan Moskowitz, Anindita Basu, Xin
He, Sebastian Pott. Single-cell genomics improves the discovery of risk
variants and genes of Atrial Fibrillation. medRxiv 2022.02.02.22270312;
doi: <https://doi.org/10.1101/2022.02.02.22270312>
