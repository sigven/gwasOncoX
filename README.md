## gwasOncoX

The goal of **gwasOncoX** is to offer an R package that simplifies the process of retrieving
variant data and track files (VCF, BED) for low-to-moderate risk variants associated with cancer, as
found in [genome-wide association studies](). The package utilizes the 
[googledrive](https://googledrive.tidyverse.org/) R package to download the pre-processed 
and documented datasets to a local cache directory provided by the user.

### Installation

`remotes::install_github('sigven/gwasOncoX')`

### Usage

The package offers (currently) five different functions, that each retrieves a specific dataset that can be of use for gene annotation purposes.

-   [`get_variants()`](https://sigven.github.io/geneGwasX/reference/get_basic.html) - retrieves basic, non-transcript-specific gene annotations. Includes tumor suppressor gene/oncogene/driver annotations from multiple resources, NCBI gene summary descriptions, as well as multiple predictions/scores when it comes to gene indispensability and loss-of-function tolerance

-   [`get_bed()`](https://sigven.github.io/geneOncoX/reference/get_gencode.html) - retrieves BED tracks ( *grch37* and *grch38* ) for variants associated with cancer

-   [`get_vcf()`](https://sigven.github.io/geneOncoX/reference/get_alias.html) - retrieves VCF files ( *grch37* and *grch38* ) for variants associated with cancer. 

### IMPORTANT NOTE

If you use the datasets provided with **gwasOncoX**, make sure you properly cite the the [NHGRI-EBI Catalog of human genome-wide association studies]()

### Contact

sigven AT ifi.uio.no
