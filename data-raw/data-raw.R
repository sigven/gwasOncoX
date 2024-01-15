suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(myvariant))

source("data-raw/gwas_utils.R")

options(timeout = 500000)

catalog_version_date <- '2023-12-20'
ebi_catalog_version_date <- '20231220'
fname_catalog_associations <- 
  file.path(
    "data-raw", 
    paste0("gwas_catalog_all_associations_", 
           ebi_catalog_version_date,".tsv"))


gwas_metadata <- 
  data.frame(
    'source' = "GWAS Catalog",
    'source_description' = "The NHGRI-EBI Catalog of human genome-wide association studies",
    'source_url' = "https://www.ebi.ac.uk/gwas/home",
    'source_citation' = "Sollis et al., Nucleic Acids Res, 2022; 36350656",	
    'source_version' = as.character(glue::glue('v{ebi_catalog_version_date}')),
    'source_abbreviation' = 'gwas_catalog',
    'source_license' = "EMBL-EBI terms of use",
    'source_license_url' = "https://www.ebi.ac.uk/about/terms-of-use"
  )

gwas_collections <- c('cancer','all')
gwas_hits_pr_rsid <- list()
gwas_hits <- list()
gwas_citations <- list()
for (c in gwas_collections) {
  gwas_hits_pr_rsid[[c]] <- data.frame()
  gwas_citations[[c]] <- data.frame()
  gwas_hits[[c]] <- NULL
}

if (!file.exists(paste0(fname_catalog_associations,".gz"))) {

  download.file(url = "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative", 
                destfile = fname_catalog_associations)
  system(paste0("gzip ", fname_catalog_associations))
}

phenotype_maps <- phenOncoX::get_aux_maps(
  cache_dir = "data-raw"
)

nuc_chromosomes_df <- data.frame("chrom" = c(as.character(seq(1:22)), "X", "Y"),
                                 stringsAsFactors = F)
all_gwas_variants <- as.data.frame(
  read.table(file = gzfile(paste0(fname_catalog_associations,".gz")), 
             sep = "\t", comment.char = "",stringsAsFactors = F, 
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
          paste(
            "(http://www.orpha.net/ORDO/)",
            "(http://purl.obolibrary.org/obo/)",
            "(http://www.ebi.ac.uk/efo/)",
            sep = "|"
          ), "")) |>
    tidyr::separate_rows(efo_id, sep = ", ") |>
    tidyr::separate_rows(snps, sep = "; ") |>
    dplyr::filter(!startsWith(efo_id,"http")) |>
    ## ignore associations to gene ontology (restrict to EFO ontology)
    #dplyr::filter(!stringr::str_detect(efo_id,"obolibrary")) |>
    dplyr::select(-c(mapped_trait_uri,intergenic,cnv)) |>
    dplyr::rename(pmid = pubmedid, rsid = snps, 
                  chromosome = chr_id, cytoband = region) |>
    dplyr::mutate(tag_snp = 'tag') |>
    dplyr::mutate(efo_id = stringr::str_replace_all(efo_id,"_",":")) |>
    dplyr::left_join(
      phenotype_maps$records$efo$efo2name, by = c("efo_id")) |>
    dplyr::mutate(gwas_hit = paste(
      rsid, strongest_snp_risk_allele, pmid, 
      tag_snp, p_value, efo_id, sep = "|"))
)

gwas_hits[['all']] <- all_gwas_variants
gwas_hits[['cancer']] <- all_gwas_variants |>
  dplyr::filter(stringr::str_detect(
    tolower(mapped_trait),
    paste0("tumor|cancer|neuroblastom|neoplasm|chemotherapy|",
           "glioma|glioblastoma|wilms|myeloma|adenocarcinoma|barrett|",
           "sarcoma|melanoma|leukaemia|leukemia|lymphom|",
           "platinum|carcinoma|hereditary|nephroblastoma"))) |> 
  dplyr::filter(
    !stringr::str_detect(
      disease_trait,
      "Select biomarker|mass index|Age-related|PCA3|Ileal|Obesity|levels")) |> 
  dplyr::filter(!stringr::str_detect(mapped_trait,"asbestos|exposure"))




gwas_hits_pr_rsid[['cancer']] <- as.data.frame(
  gwas_hits[['cancer']] |> 
    dplyr::group_by(rsid) |> 
    dplyr::summarise(gwas_hit = paste(gwas_hit, collapse = ","), 
                     .groups = "drop")
  )
gwas_hits_pr_rsid[['all']] <- as.data.frame(
  gwas_hits[['all']] |> 
    dplyr::group_by(rsid) |> 
    dplyr::summarise(gwas_hit = paste(gwas_hit, collapse = ","), 
                     .groups = "drop")
  )


version_bumped <- paste0(
  "0.",
  as.character(as.integer(substr(as.character(
    packageVersion("gwasOncoX")), 3, 3)) + 1),
  ".0")

for (c in gwas_collections) {
  
  ## GET CITATION DATA
  gwas_citations[[c]] <- get_citations_pubmed(
    pmids = unique(gwas_hits[[c]]$pmid),
    cache_pmid_fname = file.path(
      "data-raw",
      paste0(
        "citations_gwas_",c,"_current.rds"
    ))
  )
  

  ## GET DBSNP DATA (CHROM, POS, REF, ALT)
  gwas_vcf_data <- get_dbsnp_data(
    rsids = unique(gwas_hits_pr_rsid[[c]]$rsid),
    cache_dbsnp_fname = file.path(
      "data-raw",
      paste0(
        "dbsnp_gwas_",c,"_current.rds"
      ))) |>
    dplyr::left_join(
      gwas_hits_pr_rsid[[c]], by = c("rsid"),
      relationship = "many-to-many") |>
    tidyr::separate_rows(gwas_hit, sep=",") |> 
    tidyr::separate(
      gwas_hit, 
      c("rsid2","risk_allele","pmid","tag","pvalue","efo_id"),
      sep="\\|") |> 
    dplyr::filter(risk_allele != "NA") |> 
    dplyr::filter(risk_allele == ref | risk_allele == alt) |> 
    dplyr::distinct() |> 
    dplyr::mutate(
      gwas_hit = paste(
        rsid,risk_allele,pmid,"tag",pvalue,efo_id,sep="|")) |> 
    dplyr::filter(as.numeric(pvalue) <= 0.000001) |> 
    dplyr::select(-c(rsid2,risk_allele,pmid,tag,pvalue,efo_id)) |>
    dplyr::group_by(
      dplyr::across(c(-gwas_hit))
    ) |>
    dplyr::summarise(
      gwas_hit = paste(gwas_hit, collapse=","),
      .groups = "drop"
    )
  
  gwas_phenotype_records <- as.data.frame(
    gwas_hits[[c]] |>
      dplyr::left_join(
        gwas_citations[[c]],
        relationship = "many-to-many") |>
      dplyr::mutate(gwas_catalog_version = ebi_catalog_version_date) |>
      dplyr::mutate(p_value_verbose = sprintf("%2.1e",p_value)) |>
      dplyr::mutate(
        gwas_citation_exp = paste0(
          paste0('<a href=\'https://www.ebi.ac.uk/gwas/variants/',
                 rsid,'\' target=\'_blank\'>'),
          stringr::str_to_title(efo_name),"</a>, ",
          link," (association p-value = ",p_value_verbose,")")) |>
      dplyr::mutate(gwas_phenotype = stringr::str_to_title(efo_name)) |>
      dplyr::filter(!is.na(gwas_phenotype)) |>
      dplyr::select(rsid, chromosome, cytoband, dplyr::everything()) |>
      dplyr::distinct() |>
      dplyr::inner_join(
        dplyr::select(gwas_vcf_data, rsid),
        by = "rsid",
        relationship = "many-to-many"
      )
  )
  
  gwas_phenotype_data <- list()
  gwas_phenotype_data[['records']] <- gwas_phenotype_records
  gwas_phenotype_data[['metadata']] <- gwas_metadata
    
  
  
  fname_rds <- file.path(
    "data-raw", "gd_local",
    paste0("gwas_v", version_bumped, ".rds"))
  if (c == 'all') {
    fname_rds <- file.path(
      "data-raw", "gd_local",
      paste0("gwas_all_v", version_bumped, ".rds"))
  }
  
  saveRDS(gwas_phenotype_data, file = fname_rds)
  
  
  print_gwas_vcf(
    gwas_vcf_data = gwas_vcf_data, 
    subset = c, 
    output_directory = file.path(
      here::here(),
      "data-raw",
      "gd_local"
    ),
    pversion = paste0(
      "v", version_bumped))
  
  print_gwas_bed(
    gwas_vcf_data, 
    subset = c, 
    output_directory = file.path(
      here::here(),
      "data-raw",
      "gd_local"
    ),
    pversion = paste0(
      "v", version_bumped))
  
}


#googledrive::drive_auth_configure(api_key = Sys.getenv("GD_KEY"))
gd_records <- data.frame()

for (elem in c('all','cancer')) {
  
  cancer_only <- FALSE
  prefix = paste0("gwas_all_v", version_bumped)
  if (elem == "cancer") {
    prefix = paste0("gwas_v", version_bumped)
    cancer_only <- TRUE
  }
  
  for (format in c("bed", "vcf")) {
    for (build in c("grch37", "grch38")) {
      for (ftype in c("gz", "gz.tbi")) {
        (gd_rec <- googledrive::drive_upload(
          file.path("data-raw", "gd_local",
                    paste0(prefix,"_",build,".",format, ".", ftype)),
          paste0("gwasOncoX/", prefix,"_",build,".", format, ".", ftype)
        ))
        
        df_rec <- as.data.frame(gd_rec) |>
          dplyr::select(name, id) |>
          dplyr::mutate(
            datatype = format,
            assembly = build,
            cancer_only = cancer_only,
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
      cancer_only = cancer_only,
      md5Checksum = gd_rec$drive_resource[[1]]$md5Checksum)
  gd_records <- gd_records |>
    dplyr::bind_rows(
      df_rec
    )
  
  if (cancer_only == T) {
    (gd_rec <- googledrive::drive_upload(
      file.path("data-raw", "gd_local", "gwas.vcfanno.vcf_info_tags.txt"),
      paste0("gwasOncoX/gwas.vcfanno.vcf_info_tags.txt")
    ))
    
    df_rec <- as.data.frame(gd_rec) |>
      dplyr::select(name, id) |>
      dplyr::mutate(
        datatype = "vcf",
        assembly = NA,
        cancer_only = cancer_only,
        md5Checksum = gd_rec$drive_resource[[1]]$md5Checksum)
    gd_records <- gd_records |>
      dplyr::bind_rows(
        df_rec
      )
  }else{
    (gd_rec <- googledrive::drive_upload(
      file.path(
        "data-raw", 
        "gd_local", 
        "gwas_all.vcfanno.vcf_info_tags.txt"),
      paste0("gwasOncoX/gwas_all.vcfanno.vcf_info_tags.txt")
    ))
    
    df_rec <- as.data.frame(gd_rec) |>
      dplyr::select(name, id) |>
      dplyr::mutate(
        datatype = "vcf",
        assembly = NA,
        cancer_only = cancer_only,
        md5Checksum = gd_rec$drive_resource[[1]]$md5Checksum)
    gd_records <- gd_records |>
      dplyr::bind_rows(
        df_rec
      )
  }
  
  gd_records <- gd_records |>
    dplyr::mutate(date = Sys.Date(),
                  pVersion = version_bumped)
  
}

db_id_ref <- gd_records
db_id_ref$ebi_catalog_version <- ebi_catalog_version_date

usethis::use_data(db_id_ref, internal = T, overwrite = T)
