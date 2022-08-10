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
get_gwas_data <- function(cache_dir = NA,
                         overwrite = F,
                         cancer_only = T,
                         build = "grch38",
                         data_type = "vcf"){

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
  
  if(is.na(cache_dir)){
    lgr::lgr$error(paste0("Argument cache_dir = '",
                          cache_dir, "' is not defined"))
    return(0)
  }
  
  if(!dir.exists(cache_dir)){
    lgr::lgr$error(paste0("Argument cache_dir = '",
                          cache_dir, "' does not exist"))
    return(0)
  }
  
  dir_local <- file.path(
    cache_dir,
    paste0("gwasOncoX_",
           unique(gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$datatype == data_type,]$pVersion))
  )
  
  if(!dir.exists(dir_local)){
    dir.create(path = dir_local)
  }
  
  data <- NULL
  files_to_retrieve <- 
    gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$datatype == data_type & 
                            gwasOncoX:::db_id_ref$assembly == build &
                            gwasOncoX:::db_id_ref$cancer_only == cancer_only,]
  if(data_type == "vcf"){
    files_to_retrieve <- 
      gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$datatype == data_type & 
                              (is.na(gwasOncoX:::db_id_ref$assembly) |
                                 gwasOncoX:::db_id_ref$assembly == build) &
                              gwasOncoX:::db_id_ref$cancer_only == cancer_only,]
  }
  if(data_type == "rds"){
    files_to_retrieve <- 
      gwasOncoX:::db_id_ref[gwasOncoX:::db_id_ref$datatype == data_type & 
                              gwasOncoX:::db_id_ref$cancer_only == cancer_only,]
  }
  
  
  for(i in 1:nrow(files_to_retrieve)){
    fname_cache <- file.path(dir_local, 
                             stringr::str_replace(
                               files_to_retrieve[i, "name"],
                               "v[0-9]{1,}\\.[0-9]{1,}\\.[0-9]{1,}.",
                               ""
                             ))
    
    files_to_retrieve[i, "fname_cache"] <- fname_cache
    
    if(overwrite == F & file.exists(fname_cache)){
      lgr::lgr$info(paste0("File ", fname_cache, " already exists - ",
                           "set 'overwrite = TRUE' to re-download"))
      
      if(data_type == "rds"){
        data <- readRDS(file = fname_cache)
      }
    }
    else{
      fname_gd <- files_to_retrieve[i,"id"]
      #md5checksum_package <- files_to_retrieve[i,"md5Checksum"]
      
      googledrive::drive_deauth()
      
      lgr::lgr$info("Downloading remote dataset from Google Drive to 'cache_dir'")
      lgr::lgr$info(paste0("Local file name: ", fname_cache))
      dl <- googledrive::with_drive_quiet(
        googledrive::drive_download(
          fname_gd,
          path = fname_cache, overwrite = TRUE)
      )
      
      md5checksum_remote <- dl$drive_resource[[1]]$md5Checksum
      md5checksum_local <- tools::md5sum(fname_cache)
      names(md5checksum_local) <- NULL
      
      if(md5checksum_remote != md5checksum_local){
        lgr::lgr$error(paste0("md5 checksum of local file (", md5checksum_local,
                              ") is inconsistent with remote file (",
                              md5checksum_remote,")"))
      }
      
      if(data_type == "rds"){
        data <- readRDS(file = fname_cache)
      }
      
    }
  }
  
  retdata <- files_to_retrieve
  if(!is.null(data) & data_type == "rds"){
    retdata <- data
  }
  
  return(retdata)
}
