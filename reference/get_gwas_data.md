# Function that retrieves gwasOncoX data from Google Drive

Function that retrieves gwasOncoX data from Google Drive

## Usage

``` r
get_gwas_data(
  cache_dir = NA,
  force_download = F,
  cancer_only = T,
  build = "grch38",
  data_type = "vcf"
)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to force download if file exists in cache)

- cancer_only:

  retrieve data relevant for cancer phenotypes only

- build:

  genome assembly (grch37/grch38)

- data_type:

  type of data to be retrieved (bed, vcf, rds)
