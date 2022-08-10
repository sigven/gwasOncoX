#' Function that retrieves gwasOncoX data from Google Drive
#'
#' @param cache_dir Local directory for data download
#' @param overwrite Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#' @param cancer_only retrieve data relevant for cancer phenotypes only
#' @param build genome assembly (grch37/grch38)
#' @param data_type type of data to be retrieved (bed, vcf, rds)
#'
#' @keywords internal
#'
#'
get_gwas_data <- function(cache_dir = NULL,
                         overwrite = F,
                         cancer_only = T,
                         build = "grch38",
                         data_type = "vcf"){

  if(is.null(cache_dir) | !dir.exists(cache_dir)){
    lgr::lgr$error(paste0("Argument cache_dir = '",
                          cache_dir, "' is NULL or does not exist"))
  }

  
  dir_local <- file.path(
    cache_dir,
    paste0("gwasOncoX_",
           unique(gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$datatype == data_type,]$pVersion))
  )
  
  if(!dir.exists(dir_local)){
    ## use more robust try catch-ish? 
    system('mkdir ', dir_local)
  }
  
  files_to_retrieve <- 
    gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$datatype == data_type,]
  

  fname_gd <- googledrive::as_id(
    gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$name == db,]$gid)

  md5checksum_package <-
    gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$name == db,]$md5Checksum

  dat <- NULL
  if(file.exists(fname_local) & overwrite == F){
    dat <- readRDS(fname_local)
    if(!is.null(dat[['records']]) & !is.null(dat[['metadata']])){
      lgr::lgr$info(paste0(
        "Reading from cache_dir = '", cache_dir, "', argument overwrite = F"))
      lgr::lgr$info(paste0("Object 'gene_",db,"' sucessfully loaded"))
      if(db == 'gencode'){
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[['records']][['grch37']]), " records for",
          " genome build grch37"))
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[['records']][['grch38']]), " records for",
          " genome build grch38"))
      }else{
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[['records']]), " records"))
      }
    }

  }else{

    googledrive::drive_deauth()

    lgr::lgr$info("Downloading remote dataset from Google Drive to cache_dir")
    dl <- googledrive::with_drive_quiet(
      googledrive::drive_download(
        fname_gd,
        path = fname_local, overwrite = TRUE)
    )

    md5checksum_remote <- dl$drive_resource[[1]]$md5Checksum
    md5checksum_local <- tools::md5sum(fname_local)
    names(md5checksum_local) <- NULL

    if(md5checksum_remote == md5checksum_local){
      dat <- readRDS(fname_local)
      if(!is.null(dat[['records']]) & !is.null(dat[['metadata']])){

        lgr::lgr$info(paste0(
          "Reading from cache_dir = ' (", cache_dir, "'), argument overwrite = F"))
        lgr::lgr$info(paste0("Object 'gene_",db,"' sucessfully loaded"))
        lgr::lgr$info(paste0("md5 checksum is valid: ", md5checksum_remote))

        if(db == 'gencode'){
          lgr::lgr$info(paste0(
            "Retrieved n = ", nrow(dat[['records']][['grch37']]), " records for",
            " genome build grch37"))
          lgr::lgr$info(paste0(
            "Retrieved n = ", nrow(dat[['records']][['grch38']]), " records for",
            " genome build grch38"))
        }else{
          lgr::lgr$info(paste0(
            "Retrieved ", nrow(dat[['records']]), " records"))
        }
      }
    }else{
      lgr::lgr$error(paste0("md5 checksum of local file (", md5checksum_local,
                            ") is inconsistent with remote file (",
                            md5checksum_remote,")"))
    }

  }
  return(dat)
}
