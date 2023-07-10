#' @title
#' Get GWAS-associated variants
#'
#' @description
#' Downloads and returns a dataset with variants associated with
#' disease phenotypes - as discovered from genome-wide association studies.
#'
#' @param cache_dir Local directory for data download
#' @param force_overwrite Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#' @param cancer_only logical indicating retrieval of all variants versus 
#' cancer-associated variants only 
#'
#' @returns a list object with two elements
#' 1. records - A data frame with varions annotations (phenotype ID's, 
#' association p-value, risk alleles, PMIDs) for GWAS Catalog variants
#' 2. metadata - GWAS Catalog metadata
#'
#' @export
#'



get_variants <- function(cache_dir = NA,
                    force_overwrite = F,
                    cancer_only = T){
  
  vardata <- get_gwas_data(cache_dir = cache_dir,
                force_download = force_overwrite,
                data_type = "rds",
                cancer_only = cancer_only)
  return(vardata)
}