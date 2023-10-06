
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

### [Data preparation](https://xinhe-lab.github.io/mapgen/articles/data_preparation_tutorial.html)

Prepare input data: GWAS summary statistics, LD information, etc.

### [Enrichment analysis](https://xinhe-lab.github.io/mapgen/articles/enrichment_tutorial.html)

Assess the enrichment of genetic signals of a trait of interest in
functional annotations using `TORUS`.

\*Please install [TORUS](https://github.com/xqwen/torus) software
package, if you need to run enrichment analysis.

### [Fine-mapping](https://xinhe-lab.github.io/mapgen/articles/finemapping_tutorial.html)

Perform Bayesian statistical fine-mapping using `SuSiE` on
trait-associated loci, using a informative prior that favors variants
located in enriched annotations.

\*Please install [susieR](https://github.com/stephenslab/susieR)
package, if you need to run fine-mapping with GWAS summary statistics.

### [Gene mapping](https://xinhe-lab.github.io/mapgen/articles/gene_mapping_tutorial.html)

Infer causal genes at each locus based on fine-mapping result and
genomic annotations.

### [PIP partitioning by annotation categories](https://xinhe-lab.github.io/mapgen/articles/partition_pip_tutorial.html)

Partitioning finemapping PIPs by annotation categories.

### [Making track plots](https://xinhe-lab.github.io/mapgen/articles/track_plot_tutorial.html)

Making track plots of GWAS, finemapping, and annotation data using
`Gviz` package.

## Reference

> Alan Selewa\*, Kaixuan Luo\*, Michael Wasney, Linsin Smith, Xiaotong
> Sun, Chenwei Tang, Heather Eckart, Ivan Moskowitz, Anindita Basu, Xin
> He, Sebastian Pott. Single-cell genomics improves the discovery of
> risk variants and genes of atrial fibrillation. *Nat Commun.* 2023 Aug
> 17;14(1):4999. doi: 10.1038/s41467-023-40505-5. PMID: 37591828; PMCID:
> PMC10435551.
