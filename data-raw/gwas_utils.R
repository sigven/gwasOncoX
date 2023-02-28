#' A function that retrieves variant data from dbSNP (chrom,pos,ref,alt) for a a list of rsids
#'
#' @param rsids
#' @param cache_dbsnp_fname File with retrieved dbsnp data (.rds)

#' @return dbsnp_results 
#' 
get_dbsnp_data <- function(rsids, cache_dbsnp_fname = NA){
  
  rsid_df <- data.frame('rsid' = rsids)
  cache_dbsnp <- data.frame()
  
  if(!is.na(cache_dbsnp_fname)){
    if(file.exists(cache_dbsnp_fname)){
      cache_dbsnp <- as.data.frame(
        readRDS(
          file = cache_dbsnp_fname
        )
      )
      
      rsid_df <- rsid_df |>
        dplyr::anti_join(
          cache_dbsnp, by = "rsid")
    }
  }
  
  all_hits <- data.frame()
  
  if(nrow(rsid_df) == 0){
    return(cache_dbsnp)
  }else{
    all_hits <- cache_dbsnp
  }
  
  rsids <- rsid_df$rsid
  
  
  start <- 1
  stop <- min(199,length(rsids))
  while (start < length(rsids)) {
    b <- myvariant::getVariants(
      rsids[start:stop], fields = "dbsnp")
    dbsnp_results <- data.frame(
      'rsid' = b$dbsnp.rsid, 
      'chrom' = b$dbsnp.chrom, 
      'pos_start' = b$dbsnp.hg19.start,
      'pos_end' = b$dbsnp.hg19.end, 
      'ref' = b$dbsnp.ref, 
      'alt' = b$dbsnp.alt, 
      'vartype' = b$dbsnp.vartype, 
      stringsAsFactors = F) |> 
      dplyr::distinct() |>
      dplyr::filter(!is.na(chrom) & !is.na(ref) & 
                      alt != "" & ref != "" & 
                      ref != "." & alt != '.' & 
                      !is.na(alt) & 
                      !is.na(pos_start) 
                    & !is.na(pos_end))
    cat(paste0(start,"-",stop),'\n')
    all_hits <- dplyr::bind_rows(all_hits,dbsnp_results)
    
    start <- stop + 1
    stop <- min(length(rsids),start + 199)
  }
  
  all_hits <- all_hits |> 
    dplyr::inner_join(
      nuc_chromosomes_df, by = c("chrom" = "chrom"))
  
  if(!is.na(cache_dbsnp_fname)){
    saveRDS(all_hits, file = cache_dbsnp_fname)
  }
  
  return(all_hits)
  
}


#' A function that splits an array into chunks of equal size
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))



#' A function that returns a citation with first author, 
#' journal and year for a PubMed ID
#'
#' @param pmids An array of Pubmed IDs
#' @param cache_pmid_fname File with retrieved citation data (.rds)
#' @return citation PubMed citation, with first author, journal and year
#' 
get_citations_pubmed <- function(pmids, cache_pmid_fname = NA){
  
  pmid_df <- data.frame('pmid' = pmids)
  cache_citations <- data.frame()
  
  if(!is.na(cache_pmid_fname)){
    if(file.exists(cache_pmid_fname)){
      cache_citations <- as.data.frame(
        readRDS(
          file = cache_pmid_fname
        )
      )
      
      pmid_df <- pmid_df |>
        dplyr::anti_join(
          cache_citations, by = "pmid")
    }
  }
  
  all_citations <- data.frame()
  if(nrow(pmid_df) == 0){
    return(cache_citations)
  }else{
    all_citations <- cache_citations
  }
  
  pmids <- pmid_df$pmid
  
  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(pmids,ceiling(length(pmids)/100))
  j <- 0
  cat('Retrieving PubMed citations for PMID list, total length', 
      length(pmids))
  cat('\n')
  while (j < length(pmid_chunks)) {
    pmid_chunk <- pmid_chunks[[as.character(j)]]
    cat('Processing chunk ',j,' with ',length(pmid_chunk),'PMIDS')
    cat('\n')
    pmid_string <- paste(pmid_chunk,collapse = " ")
    res <- RISmed::EUtilsSummary(pmid_string, type = "esearch", 
                                 db = "pubmed", retmax = 5000)
    result <- RISmed::EUtilsGet(res)
    year <- RISmed::YearPubmed(result)
    authorlist <- RISmed::Author(result)
    pmid_list <- RISmed::PMID(result)
    i <- 1
    first_author <- c()
    while (i <= length(authorlist)) {
      first_author <- c(first_author, 
                        paste(authorlist[[i]][1,]$LastName,
                              " et al.", sep = ""))
      i <- i + 1
    }
    journal <- RISmed::ISOAbbreviation(result)
    if (length(pmid_list) == length(first_author) & 
       length(pmid_list) == length(year) & 
       length(journal) == length(pmid_list)) {
      citations <- data.frame(
        'pmid' = as.integer(pmid_list), 
        'citation' = paste(first_author, 
                           year, journal, sep = ", "), 
        stringsAsFactors = F)
      citations$link <- 
        paste0('<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/', 
               citations$pmid,'\' target=\'_blank\'>', 
               citations$citation,'</a>')
      all_citations <- dplyr::bind_rows(
        all_citations, citations)
    }else{
      cat('ERROR')
    }
    j <- j + 1
  }
  
  if(!is.na(cache_pmid_fname)){
    saveRDS(all_citations, file = cache_pmid_fname)
  }
  
  return(all_citations)
  
}

#' A function that prints GWAS data to VCF file
#'
#' @param gwas_vcf_data GWAS data
#' @param subset subset: all/cancer
#' @param output_directory output directory
#' @param pversion gwasOncoX version
#' 
print_gwas_vcf <- function(gwas_vcf_data, 
                      subset = "all", 
                      output_directory = NA,
                      pversion = "v0.2.0"){
  
  vcfanno_fname <- 
    file.path(output_directory,
              "gwas_all.vcfanno.vcf_info_tags.txt")
  
  vcf_fname_prefix <- paste(
    "gwas_all", pversion, sep = "_")
  
  if (subset == "cancer") {
    
    vcf_fname_prefix <- paste(
      "gwas", pversion, sep = "_")
    vcfanno_fname <- 
      file.path(output_directory,
                "gwas.vcfanno.vcf_info_tags.txt")
  }
  
  header_lines <- 
    c("##fileformat=VCFv4.2",
      paste0(
        "##INFO=<ID=GWAS_HIT,Number=.,Type=String,Description",
        "=\"SNP associated with disease phenotype from genome-wide ",
        "association study, format: rsid|risk_allele|pmid|tag_snp|",
        "p_value|efo_id\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

  gwas_vcf_records <- gwas_vcf_data |>
    dplyr::rename(CHROM = chrom, POS = 
                    pos_start, 
                  REF = ref, 
                  ALT = alt, 
                  INFO = gwas_hit) |>
    dplyr::filter(!is.na(REF) & !is.na(ALT) & !is.na(CHROM)) |>
    dplyr::mutate(QUAL = ".",
                  FILTER = "PASS",
                  ID = ".") |>
    dplyr::mutate(INFO = dplyr::if_else(
      !stringr::str_detect(INFO,"^GWAS_HIT="),
      paste0("GWAS_HIT=", INFO),
      as.character(INFO))) |>
    dplyr::select(
      CHROM, POS, ID, REF, 
      ALT, QUAL, FILTER, INFO) |>
    dplyr::distinct()
  
  gwas_vcf_fnames <- vcfhelpR::write_vcf_records(
    vcf_records = gwas_vcf_records,
    output_dir = output_directory,
    header_lines = header_lines,
    genome_build = "grch37",
    vcf_fname_prefix = vcf_fname_prefix,
    keep_uncompressed = TRUE
  )
  
  vcfhelpR::crossmap_vcf(
    target_genome_file = "/Users/sigven/research/DB/hg38/hg38.fa",
    direction = "hg19Tohg38",
    source_vcf = gwas_vcf_fnames[['vcf_raw']],
    target_vcf = stringr::str_replace(
      gwas_vcf_fnames[['vcf_raw']], 
      "grch37", "grch38")
  )

  write(header_lines[2],file = vcfanno_fname, sep = "\n")
}


#' A function that prints GWAS variant data to BED file
#'
#' @param gwas_vcf_data GWAS data 
#' @param subset class: all/cancer
#' @param output_directory output directory
#' @param pversion gwasOncoX version
#' 
print_gwas_bed <- function(
    gwas_vcf_data, 
    subset = "all", 
    output_directory = NA,
    pversion = "v0.6.1"){
  
  bed_fname_prefix <- paste(
    "gwas_all", pversion, sep = "_")
  
  if (subset == "cancer") {
    bed_fname_prefix <- paste(
      "gwas", pversion, sep = "_")
  }
  
  gwas_bed_records <- as.data.frame(
    gwas_vcf_data |>
      dplyr::filter(!is.na(pos_start) & !is.na(chrom)) |>
      dplyr::filter(stringr::str_detect(pos_start,"[0-9]{1,}")) |>
      dplyr::filter(nchar(pos_start) >= 4) |>
      dplyr::mutate(start = as.integer(pos_start) - 1) |>
      dplyr::rename(end = pos_start, 
                    name = gwas_hit) |>
      tidyr::separate_rows(name, sep = ",") |>
      dplyr::group_by(chrom, start, end) |>
      dplyr::summarise(name = paste(sort(unique(name)), 
                                    collapse = "@"),
                       .groups = "drop") |>
      dplyr::mutate(chrom = paste0("chr", chrom)) |>
      dplyr::mutate(chrom = stringr::str_replace(chrom,"chr",""))
  )

  bed_fnames <- vcfhelpR::write_bed_records(
    bed_records = gwas_bed_records,
    output_dir = output_directory,
    bed_fname_prefix = bed_fname_prefix,
    genome_build = "grch37",
    keep_uncompressed = T
  )
  
  vcfhelpR::crossmap_bed(
    direction = "hg19Tohg38",
    source_bed = bed_fnames[['bed_raw']],
    target_bed = stringr::str_replace(
      bed_fnames[['bed_raw']],
      "grch37","grch38"),
    remap_ratio = 1
  )

}
