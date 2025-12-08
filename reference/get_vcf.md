# Get VCF file with GWAS variants

Downloads a tabix-indexed VCF file (`.vcf.gz` and `.vcf.gz.tbi`) with
variants associated with disease phenotypes (as found in the NHGRI-EBI
GWAS catalog). Only variants where the risk-allele is properly
identified are included here. The GWAS_HIT element of the INFO column
has the following format:
rsid\|risk_allele\|pmid\|tag_snp\|p_valule\|efo_id

## Usage

``` r
get_vcf(cache_dir = NA, force_download = F, build = "grch37", cancer_only = T)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to force download if file exists in cache)

- build:

  genome assembly (grch37/grch38)

- cancer_only:

  logical indicating retrieval of all variants versus cancer-associated
  variants only

## Value

A data frame with file download information and metadata
