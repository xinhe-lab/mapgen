
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Mapgen

<!-- badges: start -->
<!-- badges: end -->

Mapgen is a multi-function software that performs the following tasks:

1.  Enrichment analysis of functional annotations for a trait of
    interest.
2.  Functionally-informed genetic fine-mapping.
3.  Gene mapping based on fine-mapping result and genomic annotations.

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

## Tutorials

### 1. [Data preparation](https://xinhe-lab.github.io/mapgen/articles/data_preparation_tutorial.html)

Prepare input data: GWAS summary statistics, LD reference panel, etc.

### 2. [Enrichment analysis](https://xinhe-lab.github.io/mapgen/articles/enrichment_finemapping_tutorial.html)

Assess the enrichment of genetic signals of a trait of interest in
functional annotations using `TORUS`.

\*Please install [TORUS](https://github.com/xqwen/torus) software
package, if you need to run enrichment analysis.

### 3. [Fine-mapping](https://xinhe-lab.github.io/mapgen/articles/enrichment_finemapping_tutorial.html)

Perform Bayesian statistical fine-mapping using `SuSiE` on
trait-associated loci, using a informative prior that favors variants
located in enriched annotations.

\*Please install [susieR](https://github.com/stephenslab/susieR)
package, if you need to run fine-mapping with GWAS summary statistics.

### 4. [Gene mapping](https://xinhe-lab.github.io/mapgen/articles/gene_mapping_tutorial.html)

Infer causal genes at each locus based on fine-mapping result and
genomic annotations, including gene annotations, chromatin loops, etc.

## Reference

Alan Selewa\*, Kaixuan Luo\*, Michael Wasney, Linsin Smith, Xiaotong
Sun, Chenwei Tang, Heather Eckart, Ivan Moskowitz, Anindita Basu, Xin
He, Sebastian Pott. Single-cell genomics improves the discovery of risk
variants and genes of Atrial Fibrillation. medRxiv 2022.02.02.22270312;
doi: <https://doi.org/10.1101/2022.02.02.22270312>
