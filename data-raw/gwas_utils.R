#' A function that retrieves variant data from dbSNP (chrom,pos,ref,alt) for a a list of rsids
#'
#' @param rsids
#' @return dbsnp_results 
#' 
get_dbsnp_data <- function(rsids){
  
  all_hits <- data.frame()
  start <- 1
  stop <- min(199,length(rsids))
  while (start < length(rsids)) {
    b <- myvariant::getVariants(rsids[start:stop], fields = "dbsnp")
    dbsnp_results <- data.frame('rsid' = b$dbsnp.rsid, 
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
    dplyr::inner_join(nuc_chromosomes_df, by = c("chrom" = "chrom"))
  
  return(all_hits)
  
}


#' A function that splits an array into chunks of equal size
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))



#' A function that returns a citation with first author, 
#' journal and year for a PubMed ID
#'
#' @param pmid An array of Pubmed IDs
#' @return citation PubMed citation, with first author, journal and year
#' 
get_citations_pubmed <- function(pmid){
  
  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(pmid,ceiling(length(pmid)/100))
  j <- 0
  all_citations <- data.frame()
  cat('Retrieving PubMed citations for PMID list, total length', length(pmid))
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
                        paste(authorlist[[i]][1,]$LastName," et al.", sep = ""))
      i <- i + 1
    }
    journal <- RISmed::ISOAbbreviation(result)
    if (length(pmid_list) == length(first_author) & 
       length(pmid_list) == length(year) & 
       length(journal) == length(pmid_list)) {
      citations <- data.frame('pmid' = as.integer(pmid_list), 
                              'citation' = paste(first_author, 
                                                 year, journal, sep = ", "), 
                              stringsAsFactors = F)
      citations$link <- 
        paste0('<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/', 
               citations$pmid,'\' target=\'_blank\'>', 
               citations$citation,'</a>')
      all_citations <- dplyr::bind_rows(all_citations, citations)
    }else{
      cat('ERROR')
    }
    j <- j + 1
  }
  
  return(all_citations)
  
}

#' A function that prints GWAS data to VCF file
#'
#' @param gwas_vcf_data GWAS data
#' @param cl class: all/cancer
#' @param pversion gwasOncoX version
#' 
print_vcf <- function(gwas_vcf_data, cl = "all", pversion = "v0.2.0"){
  
  vcf_fname <- 
    file.path("data-raw", 
              "gd_local",
              paste0("gwas_all.",pversion,".grch37.vcf"))
  vcf_fname_grch38 <- 
    file.path("data-raw", 
              "gd_local",
              paste0("gwas_all.",pversion,".grch38.vcf"))
  vcf_content_fname <- 
    file.path("data-raw", 
              "gwas_all_vcfcontent.vcf")
  vcfanno_fname <- 
    file.path("data-raw", 
              "gd_local",
              "gwas_all.vcfanno.vcf_info_tags.txt")
  
  if (cl == "cancer") {
    vcf_fname <- 
      file.path("data-raw", 
                "gd_local",
                paste0("gwas.", 
                       pversion, ".grch37.vcf"))
    vcf_fname_grch38 <- 
      file.path("data-raw", 
                "gd_local",
                paste0("gwas.", 
                       pversion, ".grch38.vcf"))
    vcf_content_fname <- 
      file.path("data-raw", 
                "gwas_vcfcontent.vcf")
    vcfanno_fname <- 
      file.path("data-raw", 
                "gd_local",
                "gwas.vcfanno.vcf_info_tags.txt")
  }
  
  header_lines <- 
    c("##fileformat=VCFv4.2","##assembly=grch37", 
      "##INFO=<ID=GWAS_HIT,Number=.,Type=String,Description=\"SNP associated with disease phenotype from genome-wide association study, format: rsid|risk_allele|pmid|tag_snp|p_value|efo_id\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  write(header_lines, file = vcf_fname, sep = "\n")
  
  gwas_vcf <- gwas_vcf_data |>
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
    dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) |>
    dplyr::distinct()
  

  write.table(gwas_vcf, file = vcf_content_fname,
              sep = "\t", col.names = F, quote = F, 
              row.names = F)
  
  system(paste0(
    "cat ",vcf_content_fname,
    " | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5 >> ",
    vcf_fname))
  system(paste0("cat ",vcf_content_fname, 
                " | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ",
                vcf_fname))
  system(paste0("bgzip -c ", vcf_fname," > ",vcf_fname, ".gz"))
  system(paste0("tabix -p vcf ", vcf_fname, ".gz"))
  
  write(header_lines[3],file = vcfanno_fname, sep = "\n")
  
  crossmapr::crossmap_vcf(
    target_genome_file = '/Users/sigven/research/DB/hg38/hg38.fa', 
    direction = 'hg19Tohg38', 
    source_vcf = vcf_fname, 
    target_vcf = vcf_fname_grch38)
  system(paste0("rm -f ", vcf_content_fname))
  system(paste0("rm -f ", vcf_fname))
  system(paste0("rm -f ", vcf_fname_grch38))
}


#' A function that prints GWAS variant data to BED file
#'
#' @param gwas_vcf_data GWAS data 
#' @param cl class: all/cancer
#' @param pversion gwasOncoX version
#' 
print_bed <- function(gwas_vcf_data, cl = "all", pversion = "v0.2.0"){
  
  bed_fname <- 
    file.path("data-raw", 
              "gd_local",
              paste0("gwas_all.",pversion,".grch37.bed"))
  bed_fname_grch38 <- 
    file.path("data-raw", 
              "gd_local",
              paste0("gwas_all.",pversion,".grch38.bed"))
  bed_content_fname <- 
    file.path("data-raw", 
              "gwas_all_bedcontent.bed")
  
  if (cl == "cancer") {
    bed_fname <- 
      file.path("data-raw", 
                "gd_local",
                paste0("gwas.", 
                       pversion, ".grch37.bed"))
    bed_fname_grch38 <- 
      file.path("data-raw", 
                "gd_local",
                paste0("gwas.", 
                       pversion, ".grch38.bed"))
    bed_content_fname <- 
      file.path("data-raw", 
                "gwas_bedcontent.bed")
  }
  
  gwas_bed <- as.data.frame(
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
      sort_bed_regions() |>
      dplyr::mutate(chrom = stringr::str_replace(chrom,"chr",""))
  )

  options(scipen = 999)
  write.table(gwas_bed, file = bed_fname, 
              sep = "\t",col.names = F, quote = F, row.names = F)
  
  system(paste0("bgzip -c ",bed_fname," > ",bed_fname,".gz"))
  system(paste0("tabix -p bed ",bed_fname,".gz"))
  
  tmp <- as.data.frame(
    readr::read_tsv(bed_fname, col_names = F, show_col_types = F))
  tmp$X1 <- paste0('chr',tmp$X1)
  write.table(tmp, file = paste0("data-raw/gwas_all.cm.bed"), 
              col.names = F, row.names = F, quote = F, sep = "\t")
  
  crossmapr::crossmap_bed(direction = 'hg19Tohg38', 
                          source_bed = paste0("data-raw/gwas_all.cm.bed"), 
                          target_bed = bed_fname_grch38)
  system("rm -f data-raw/gwas_all.cm.bed")
  system(paste0("rm -f ", bed_fname))
  system(paste0("rm -f ", bed_fname_grch38))
}

sort_bed_regions <- function(unsorted_regions){
  sorted_regions <- NULL
  assertable::assert_colnames(
    unsorted_regions,c('start','end'),only_colnames = F, quiet = T)
  if ("chrom" %in% colnames(unsorted_regions) & 
     "start" %in% colnames(unsorted_regions) & 
     "end" %in% colnames(unsorted_regions)) {
    chrOrder <- paste0('chr',c(as.character(c(1:22)), "X", "Y", "M"))
    unsorted_regions$chrom <- factor(unsorted_regions$chrom, 
                                     levels = chrOrder)
    unsorted_regions <- unsorted_regions[order(unsorted_regions$chrom), ]
    
    sorted_regions <- data.frame()
    for (chrom in chrOrder) {
      if (nrow(unsorted_regions[unsorted_regions$chrom == chrom,]) > 0) {
        chrom_regions <- unsorted_regions[unsorted_regions$chrom == chrom,]
        chrom_regions_sorted <- 
          chrom_regions[with(chrom_regions, order(start, end)),]
        sorted_regions <- 
          dplyr::bind_rows(sorted_regions, chrom_regions_sorted)
      }
    }
    
    sorted_regions$start <- as.integer(sorted_regions$start)
    sorted_regions$end <- as.integer(sorted_regions$end)
  }
  return(sorted_regions)
  
}

