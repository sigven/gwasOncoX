#' @title
#' Get BED file with GWAS variants
#'
#' @description
#' Downloads a tabix-indexed BED file (`.bed.gz` and `.bed.gz.tbi`) with 
#' low to modest risk variants  associated with disease phenotypes (as found in 
#' the NHGRI-EBI GWAS catalog). Only variants where the risk-allele is properly 
#' identified are included here. The name of each record in the BED file has 
#' the following format: rsid|risk_allele|pmid|tag_snp|p_valule|efo_id 
#' (multiple instances separated by the `at` sign)
#'
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
#' (set to TRUE to force download if file exists in cache)
#' @param build genome assembly (grch37/grch38)
#' @param cancer_only logical indicating retrieval of all variants versus 
#' cancer-associated variants only 
#'
#' @returns
#' A data frame with file download information and metadata
#'
#' @export
#'

get_bed <- function(cache_dir = NA,
                    force_download = F,
                    build = "grch37",
                    cancer_only = T){
  
  stopifnot(build == "grch38" | build == "grch37")

  file_metadata <- get_gwas_data(cache_dir = cache_dir,
                force_download = force_download,
                data_type = "bed",
                cancer_only = cancer_only,
                build = build)
  
  return(file_metadata)
  
}

