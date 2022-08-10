suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(myvariant))

source("data-raw/gwas_utils.R")

options(timeout=50000)

catalog_version_date <- '2022-07-30'
ebi_catalog_version_date <- '20220730'

gwas_collections <- c('cancer','all')
gwas_hits_pr_rsid <- list()
gwas_hits <- list()
gwas_citations <- list()
for(c in gwas_collections){
  gwas_hits_pr_rsid[[c]] <- data.frame()
  gwas_citations[[c]] <- data.frame()
  gwas_hits[[c]] <- NULL
}

download.file(url = "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative", 
              destfile = file.path(
                "data-raw", 
                paste0("gwas_catalog_all_associations_",ebi_catalog_version_date,".tsv")))
system(paste0("gzip data-raw/gwas_catalog_all_associations_",ebi_catalog_version_date,".tsv"))

nuc_chromosomes_df <- data.frame("chrom" = c(as.character(seq(1:22)), "X", "Y"),
                                 stringsAsFactors = F)
all_gwas_variants <- as.data.frame(
  read.table(file = gzfile(paste0("data-raw/gwas_catalog_all_associations_", 
                         ebi_catalog_version_date,".tsv.gz")), 
             sep="\t", comment.char = "",stringsAsFactors = F, 
             header = T,quote = "", fill = T) |>
    janitor::clean_names() |>
    dplyr::select(pubmedid, disease_trait, region, 
                  chr_id,snps, strongest_snp_risk_allele, 
                  intergenic, p_value, cnv, 
                  mapped_trait, mapped_trait_uri) |> 
    dplyr::distinct() |>
    dplyr::inner_join(nuc_chromosomes_df, by = c("chr_id" = "chrom")) |>
    dplyr::filter(startsWith(snps,"rs")) |>
    dplyr::mutate(
      strongest_snp_risk_allele = 
        dplyr::if_else(stringr::str_detect(strongest_snp_risk_allele,"-") & 
                         startsWith(strongest_snp_risk_allele, snps),
                       stringr::str_split_fixed(strongest_snp_risk_allele, 
                                                "-",n = 2)[,2],
                       as.character(NA))) |>
    dplyr::mutate(
      strongest_snp_risk_allele = 
        dplyr::if_else(stringr::str_detect(strongest_snp_risk_allele, 
                                           "^(A|C|G|T){1,}$"),
                       strongest_snp_risk_allele,
                       as.character(NA))) |>
    dplyr::filter(!stringr::str_detect(snps," x ")) |>
    dplyr::mutate(
      efo_id = 
        stringr::str_replace_all(
          mapped_trait_uri,
          "(http://www.orpha.net/ORDO/)|(http://purl.obolibrary.org/obo/)|(http://www.ebi.ac.uk/efo/)","")) |>
    tidyr::separate_rows(efo_id,sep=", ") |>
    tidyr::separate_rows(snps,sep="; ") |>
    dplyr::filter(!startsWith(efo_id,"http")) |>
    ## ignore associations to gene ontology (restrict to EFO ontology)
    #dplyr::filter(!stringr::str_detect(efo_id,"obolibrary")) |>
    dplyr::select(-c(mapped_trait_uri,intergenic,cnv)) |>
    dplyr::rename(pmid = pubmedid, rsid = snps, 
                  chromosome = chr_id, cytoband = region) |>
    dplyr::mutate(tag_snp = 'tag') |>
    dplyr::mutate(efo_id = stringr::str_replace_all(efo_id,"_",":")) |>
    dplyr::left_join(oncoPhenoMap::auxiliary_maps$efo$efo2name, by = c("efo_id")) |>
    dplyr::mutate(GWAS_HIT = paste(
      rsid, strongest_snp_risk_allele, pmid, 
      tag_snp, p_value, efo_id, sep="|"))
)

gwas_hits[['all']] <- all_gwas_variants
gwas_hits[['cancer']] <- all_gwas_variants |>
  dplyr::filter(stringr::str_detect(
    tolower(mapped_trait),"tumor|cancer|neuroblastom|neoplasm|chemotherapy|glioma|glioblastoma|wilms|myeloma|adenocarcinoma|barrett|sarcoma|melanoma|leukaemia|leukemia|lymphom|platinum|carcinoma|hereditary")) |> 
  dplyr::filter(
    !stringr::str_detect(
      disease_trait,
      "Select biomarker|mass index|Age-related|PCA3|Ileal|Obesity|levels")) |> 
  dplyr::filter(!stringr::str_detect(mapped_trait,"asbestos|exposure"))
  
gwas_hits_pr_rsid[['cancer']] <- as.data.frame(
  gwas_hits[['cancer']] |> 
    dplyr::group_by(rsid) |> 
    dplyr::summarise(GWAS_HIT = paste(GWAS_HIT, collapse=","), .groups = "drop")
  )
gwas_hits_pr_rsid[['all']] <- as.data.frame(
  gwas_hits[['all']] |> 
    dplyr::group_by(rsid) |> 
    dplyr::summarise(GWAS_HIT = paste(GWAS_HIT, collapse=","), .groups = "drop")
  )


version_minor_bumped <- paste0(
  "0.",
  as.character(as.integer(substr(as.character(packageVersion("gwasOncoX")),3,3)) + 1),
  ".0")

for(c in gwas_collections){
  
  ## GET CITATION DATA
  gwas_citations[[c]] <- get_citations_pubmed(unique(gwas_hits[[c]]$pmid))
  
  ## GET DBSNP DATA (CHROM, POS, REF, ALT)
  vcf_data <- get_dbsnp_data(unique(gwas_hits_pr_rsid[[c]]$rsid)) |>
    dplyr::left_join(gwas_hits_pr_rsid[[c]],by=c("rsid")) |>
    dplyr::mutate(GWAS_HIT = paste0('GWAS_HIT=',GWAS_HIT))
  
  gwas_data <- as.data.frame(gwas_hits[[c]] |>
    dplyr::full_join(dplyr::select(vcf_data, rsid, chrom, 
                                   pos_start, pos_end, 
                                   ref, alt, vartype)) |>
    dplyr::mutate(var_id = paste(chrom, pos_start, ref, alt, sep = "_")) |>
    dplyr::left_join(gwas_citations[[c]]) |>
    dplyr::mutate(gwas_catalog_version = ebi_catalog_version_date) |>
    dplyr::select(-chromosome) |>
    dplyr::select(var_id, chrom, pos_start, ref, alt, vartype,
                  rsid, cytoband, strongest_snp_risk_allele,
                  mapped_trait, efo_id, efo_name, primary_site,
                  dplyr::everything()) 
  )
  
  
  fname_rds <- file.path("data-raw", "gd_local",
                         paste0("gwas.v", version_minor_bumped,".rds"))
  if(c == 'all'){
    fname_rds <- file.path("data-raw", "gd_local",
                           paste0("gwas_all.v", version_minor_bumped,".rds"))
  }
  saveRDS(gwas_data, file = fname_rds)
  
  
  print_vcf(gwas_data, cl = c, pversion = paste0("v",version_minor_bumped))
  print_bed(gwas_data, cl = c, pversion = paste0("v",version_minor_bumped))
  
}


googledrive::drive_auth_configure(api_key = Sys.getenv("GD_KEY"))
gd_records <- data.frame()

for(elem in c('all','cancer')){
  
  prefix = paste0("gwas_all.v",version_minor_bumped)
  if(elem == "cancer"){
    prefix = paste0("gwas.v",version_minor_bumped)
  }
  
  for(format in c('bed','vcf')){
    for(build in c('grch37','grch38')){
      for(ftype in c('gz','gz.tbi')){
        (gd_rec <- googledrive::drive_upload(
          file.path("data-raw", "gd_local",
                    paste0(prefix,".",build,".",format,".",ftype)),
          paste0("gwasOncoX/", prefix,".",build,".", format, ".",ftype)
        ))
        
        df_rec <- as.data.frame(gd_rec) |>
          dplyr::select(name, id) |>
          dplyr::mutate(
            datatype = format,
            assembly = build,
            md5Checksum = gd_rec$drive_resource[[1]]$md5Checksum)
        
          
        gd_records <- gd_records |>
          dplyr::bind_rows(
            df_rec
          )
      }
    }
  }
  
  (gd_rec <- googledrive::drive_upload(
    file.path("data-raw", "gd_local",
              paste0(prefix,".rds")),
    paste0("gwasOncoX/", prefix,".rds")
  ))
  
  df_rec <- as.data.frame(gd_rec) |>
    dplyr::select(name, id) |>
    dplyr::mutate(
      datatype = "rds",
      assembly = NA,
      md5Checksum = gd_rec$drive_resource[[1]]$md5Checksum)
  gd_records <- gd_records |>
    dplyr::bind_rows(
      df_rec
    )
  
  (gd_rec <- googledrive::drive_upload(
    file.path("data-raw", "gd_local","gwas.vcfanno.vcf_info_tags.txt"),
    paste0("gwasOncoX/gwas.vcfanno.vcf_info_tags.txt")
  ))
  
  df_rec <- as.data.frame(gd_rec) |>
    dplyr::select(name, id) |>
    dplyr::mutate(
      datatype = "vcf",
      assembly = NA,
      md5Checksum = gd_rec$drive_resource[[1]]$md5Checksum)
  gd_records <- gd_records |>
    dplyr::bind_rows(
      df_rec
    )
  
  gd_records <- gd_records |>
    dplyr::mutate(date = Sys.Date(),
                  pVersion = version_minor_bumped)
  
}

gd_records <- 

db_id_ref <- dplyr::bind_rows(
  dplyr::select(as.data.frame(gd_records$basic), name, id),
  dplyr::select(as.data.frame(gd_records$alias), name, id),
  dplyr::select(as.data.frame(gd_records$panels), name, id),
  dplyr::select(as.data.frame(gd_records$gencode), name, id),
  dplyr::select(as.data.frame(gd_records$predisposition), name, id)) |>
  dplyr::rename(gid = id,
                filename = name) |>
  dplyr::mutate(name = stringr::str_replace(
    stringr::str_replace(filename,"_v\\S+$",""),
    "gene_","")) |>
  dplyr::mutate(date = Sys.Date(),
                pVersion = version_minor_bumped)
db_id_ref$md5Checksum <- NA
for(elem in c('basic','predisposition','panels','alias','gencode')){
  
  db_id_ref[db_id_ref$name == elem,]$md5Checksum <-
    gd_records[[elem]]$drive_resource[[1]]$md5Checksum
}

usethis::use_data(db_id_ref, internal = T, overwrite = T)
