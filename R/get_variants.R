#' @title
#' Get GWAS-associated variants
#'
#' @description
#' Downloads and returns a dataset with variants associated with
#' disease phenotypes - as discovered from genome-wide association studies.
#'
#' @param cache_dir Local directory for data download
#' @param overwrite Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#' @param cancer_only logical indicating retrieval of all variants versus 
#' cancer-associated variants only 
#'
#' @returns a data frame with varions annotations (phenotype ID's, 
#' association p-value, risk alleles, PMIDs) for GWAS variants
#'
#' @export
#'



get_variants <- function(cache_dir = NA,
                    overwrite = F,
                    cancer_only = T){
  
  vardata <- get_gwas_data(cache_dir = cache_dir,
                overwrite = overwrite,
                data_type = "rds",
                cancer_only = cancer_only)
  return(vardata)
}