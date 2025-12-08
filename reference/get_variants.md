# Get GWAS-associated variants

Downloads and returns a dataset with low to modest risk variants
associated with disease phenotypes - as discovered from genome-wide
association studies (NHGRI-EBI GWAS catalog)

## Usage

``` r
get_variants(cache_dir = NA, force_overwrite = F, cancer_only = T)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_overwrite:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

- cancer_only:

  logical indicating retrieval of all variants versus cancer-associated
  variants only

## Value

a list object with two elements

1.  records - A data frame with variant annotations (phenotype ID,
    association p-value, risk allele, PMID etc.) for GWAS Catalog
    variants

2.  metadata - GWAS Catalog metadata
