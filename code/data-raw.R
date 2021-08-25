suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(myvariant))
suppressPackageStartupMessages(library(janitor))

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - ", level, " - ", ..., "\n", collapse = "")
}

log4r_logger <- log4r::logger(threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))


ebi_catalog_version_date <- '20210816'
datestamp <- ebi_catalog_version_date

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
              destfile = paste0("data-raw/gwas_catalog_all_associations_",ebi_catalog_version_date,".tsv"))

all_gwas_variants <- as.data.frame(
  read.table(file=paste0("data-raw/gwas_catalog_all_associations_", 
                         ebi_catalog_version_date,".tsv"), 
             sep="\t", comment.char = "",stringsAsFactors = F, 
             header = T,quote = "", fill = T) %>%
  janitor::clean_names() %>%
  dplyr::select(pubmedid, disease_trait, region, 
                chr_id,snps, strongest_snp_risk_allele, 
                intergenic, p_value, cnv, 
                mapped_trait, mapped_trait_uri) %>% 
  dplyr::distinct() %>%
  pcgrr::get_ordinary_chromosomes(chrom_var = "chr_id") %>%
  dplyr::filter(startsWith(snps,"rs")) %>%
  dplyr::mutate(
    strongest_snp_risk_allele = 
      dplyr::if_else(stringr::str_detect(strongest_snp_risk_allele,"-") & 
                       startsWith(strongest_snp_risk_allele, snps),
                     stringr::str_split_fixed(strongest_snp_risk_allele, 
                                              "-",n = 2)[,2],
                     as.character(NA))) %>%
  dplyr::mutate(
    strongest_snp_risk_allele = 
      dplyr::if_else(stringr::str_detect(strongest_snp_risk_allele, 
                                         "^(A|C|G|T){1,}$"),
                     strongest_snp_risk_allele,
                     as.character(NA))) %>%
  dplyr::filter(!stringr::str_detect(snps," x ")) %>%
  dplyr::mutate(
    efo_id = 
      stringr::str_replace_all(
        mapped_trait_uri,
        "(http://www.orpha.net/ORDO/)|(http://purl.obolibrary.org/obo/)|(http://www.ebi.ac.uk/efo/)","")) %>%
  tidyr::separate_rows(efo_id,sep=", ") %>%
  tidyr::separate_rows(snps,sep="; ") %>%
  dplyr::filter(!startsWith(efo_id,"http")) %>%
  ## ignore associations to gene ontology (restrict to EFO ontology)
  #dplyr::filter(!stringr::str_detect(efo_id,"obolibrary")) %>%
  dplyr::select(-c(mapped_trait_uri,intergenic,cnv)) %>%
  dplyr::rename(pmid = pubmedid, rsid = snps, 
                chromosome = chr_id, cytoband = region) %>%
  dplyr::mutate(tag_snp = 'tag') %>%
  dplyr::mutate(efo_id = stringr::str_replace_all(efo_id,"_",":")) %>%
  dplyr::left_join(oncoPhenoMap::efo2name, by = c("efo_id")) %>%
  dplyr::mutate(GWAS_HIT = paste(rsid, strongest_snp_risk_allele, pmid, 
                                 tag_snp,p_value,efo_id,sep="|"))
)

gwas_hits[['all']] <- all_gwas_variants
gwas_hits[['cancer']] <- all_gwas_variants %>%
  dplyr::filter(stringr::str_detect(
    tolower(mapped_trait),"tumor|cancer|neuroblastom|neoplasm|chemotherapy|glioma|glioblastoma|wilms|myeloma|adenocarcinoma|barrett|sarcoma|melanoma|leukaemia|leukemia|lymphom|platinum|carcinoma|hereditary")) %>% 
  dplyr::filter(
    !stringr::str_detect(
      disease_trait,
      "Select biomarker|mass index|Age-related|PCA3|Ileal|Obesity|levels")) %>% 
  dplyr::filter(!stringr::str_detect(mapped_trait,"asbestos|exposure"))
  
gwas_hits_pr_rsid[['cancer']] <- as.data.frame(
  gwas_hits[['cancer']] %>% 
    dplyr::group_by(rsid) %>% 
    dplyr::summarise(GWAS_HIT = paste(GWAS_HIT, collapse=","), .groups = "drop")
  )
gwas_hits_pr_rsid[['all']] <- as.data.frame(
  gwas_hits[['all']] %>% 
    dplyr::group_by(rsid) %>% 
    dplyr::summarise(GWAS_HIT = paste(GWAS_HIT, collapse=","), .groups = "drop")
  )



#' A function that retrieves variant data from dbSNP (chrom,pos,ref,alt) for a a list of rsids
#'
#' @param rsids
#' @return dbsnp_results 
#' @export
get_dbsnp_data <- function(rsids){
  
  all_hits <- data.frame()
  start <- 1
  stop <- min(199,length(rsids))
  while(start < length(rsids)){
    b <- myvariant::getVariants(rsids[start:stop],fields="dbsnp")
    dbsnp_results <- data.frame('rsid' = b$dbsnp.rsid, 
                                'chrom' = b$dbsnp.chrom, 
                                'pos_start' = b$dbsnp.hg19.start,
                                'pos_end' = b$dbsnp.hg19.end, 
                                'ref' = b$dbsnp.ref, 
                                'alt' = b$dbsnp.alt, 
                                'vartype' = b$dbsnp.vartype, 
                                stringsAsFactors = F) %>% 
      dplyr::distinct() %>%
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
  all_hits <- pcgrr::get_ordinary_chromosomes(all_hits, chrom_var = "chrom")
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
  while(j < length(pmid_chunks)){
    pmid_chunk <- pmid_chunks[[as.character(j)]]
    cat('Processing chunk ',j,' with ',length(pmid_chunk),'PMIDS')
    cat('\n')
    pmid_string <- paste(pmid_chunk,collapse = " ")
    res <- RISmed::EUtilsSummary(pmid_string, type="esearch", 
                                 db="pubmed", retmax = 5000)
    result <- RISmed::EUtilsGet(res)
    year <- RISmed::YearPubmed(result)
    authorlist <- RISmed::Author(result)
    pmid_list <- RISmed::PMID(result)
    i <- 1
    first_author <- c()
    while(i <= length(authorlist)){
      first_author <- c(first_author, 
                        paste(authorlist[[i]][1,]$LastName," et al.",sep=""))
      i <- i + 1
    }
    journal <- RISmed::ISOAbbreviation(result)
    if(length(pmid_list) == length(first_author) & 
       length(pmid_list) == length(year) & 
       length(journal) == length(pmid_list)){
      citations <- data.frame('pmid' = as.integer(pmid_list), 
                              'citation' = paste(first_author, 
                                                 year, journal,sep=", "), 
                              stringsAsFactors = F)
      citations$link <- 
        paste0('<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/', 
               citations$pmid,'\' target=\'_blank\'>', 
               citations$citation,'</a>')
      all_citations <- dplyr::bind_rows(all_citations, citations)
    }else{
      cat('BALLE')
    }
    j <- j + 1
  }

  return(all_citations)
  
}

#' A function that prints VCF data to VCF file
#'
#' @param vcf_data GWAS vcf data (chrom,pos,ref,alt,info)
#' @param class all/cancer
#' @export
print_vcf <- function(vcf_data, class = "all"){
  
  vcf_fname <- paste0("data/gwas_all.",datestamp,".vcf")
  vcf_fname_grch38 <- "data/gwas_all.grch38.vcf"
  vcf_content_fname <- paste0("data/gwas_all_vcfcontent.", 
                              datestamp,".tsv")
  vcfanno_fname <- paste0("data/gwas_all.", 
                          datestamp,".vcfanno.vcf_info_tags.txt")
  vcfanno_fname_link <- "gwas_all.grch37.vcfanno.vcf_info_tags.txt"
  vcfanno_fname_link_grch38 <- "gwas_all.grch38.vcfanno.vcf_info_tags.txt"
  vcf_fname_link <- "gwas_all.grch37.vcf"
  vcf_fname_link_grch38 <- "gwas_all.grch38.vcf"
  vcf_fname_link_grch38_gz <- "gwas_all.grch38.vcf.gz"
  vcf_fname_link_grch38_gz_tbi <- "gwas_all.grch38.vcf.gz.tbi"
  vcf_fname_link_gz <- "gwas_all.grch37.vcf.gz"
  vcf_fname_link_gz_tbi <- "gwas_all.grch37.vcf.gz.tbi"
  if(class == "cancer"){
    vcf_fname <- paste0("data/gwas.",datestamp,".vcf")
    vcf_fname_grch38 <- "data/gwas.grch38.vcf"
    vcf_content_fname <- paste0("data/gwas_vcfcontent.", datestamp,".tsv")
    vcfanno_fname <- paste0("data/gwas.",datestamp,".vcfanno.vcf_info_tags.txt")
    vcfanno_fname_link <- "gwas.grch37.vcfanno.vcf_info_tags.txt"
    vcfanno_fname_link_grch38 <- "gwas.grch38.vcfanno.vcf_info_tags.txt"
    vcf_fname_link <- "gwas.grch37.vcf"
    vcf_fname_link_grch38_gz <- "gwas.grch38.vcf.gz"
    vcf_fname_link_grch38_gz_tbi <- "gwas.grch38.vcf.gz.tbi"
    vcf_fname_link_grch38 <- "gwas.grch38.vcf"
    vcf_fname_link_gz <- "gwas.grch37.vcf.gz"
    vcf_fname_link_gz_tbi <- "gwas.grch37.vcf.gz.tbi"
  }
  
  header_lines <- 
    c("##fileformat=VCFv4.2","##assembly=grch37", 
      paste0("##SOURCE_GWAS=",datestamp), 
      "##INFO=<ID=GWAS_HIT,Number=.,Type=String,Description=\"SNP associated with disease phenotype from genome-wide association study, format: rsid|risk_allele|pmid|tag_snp|p_valule|efo_id\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  write(header_lines,file=vcf_fname,sep="\n")
  
  gwas_vcf <- vcf_data
  gwas_vcf$QUAL <- '.'
  gwas_vcf$FILTER <- 'PASS'
  gwas_vcf$ID <- '.'
  gwas_vcf <-  dplyr::rename(gwas_vcf, CHROM = chrom, POS = 
                               pos_start, REF = ref, ALT = alt, 
                             INFO = GWAS_HIT)
  gwas_vcf <- gwas_vcf[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")]
  
  write.table(gwas_vcf, file=vcf_content_fname,sep="\t",col.names = F,quote=F, row.names = F)
  
  system(paste0("cat ",vcf_content_fname," | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5 >> ",vcf_fname))
  system(paste0("cat ",vcf_content_fname," | egrep \"^[XYM]\" | sort -k1,1 -k2,2n -k4,4 -k5,5 >> ",vcf_fname))
  system(paste0("bgzip -c ",vcf_fname," > ",vcf_fname,".gz"))
  system(paste0("tabix -p vcf ",vcf_fname,".gz"))
  system(paste0("ln -F -s ",vcf_fname," ",vcf_fname_link))
  system(paste0("ln -F -s ",vcf_fname,".gz ",vcf_fname_link_gz))
  system(paste0("ln -F -s ",vcf_fname,".gz.tbi ",vcf_fname_link_gz_tbi))
  
  write(header_lines[3:4],file=vcfanno_fname,sep="\n")
  system(paste0("ln -F -s ",vcfanno_fname," ",vcfanno_fname_link))
  
  crossmapr::crossmap_vcf(target_genome_file = '/Users/sigven/research/DB/hg38/hg38.fa', 
                          direction = 'hg19Tohg38', source_vcf = vcf_fname, 
                          target_vcf = vcf_fname_grch38)
  system(paste0("ln -F -s ",vcf_fname_grch38," ",vcf_fname_link_grch38))
  system(paste0("ln -F -s ",vcf_fname_grch38,".gz ",vcf_fname_link_grch38_gz))
  system(paste0("ln -F -s ",vcf_fname_grch38,".gz.tbi ",vcf_fname_link_grch38_gz_tbi))
  system(paste0("ln -F -s ",vcfanno_fname," ", vcfanno_fname_link_grch38))
}


for(c in gwas_collections){
  
  ## GET CITATION DATA
  gwas_citations[[c]] <- get_citations_pubmed(unique(gwas_hits[[c]]$pmid))
  
  ## GET DBSNP DATA (CHROM, POS, REF, ALT)
  vcf_data <- get_dbsnp_data(unique(gwas_hits_pr_rsid[[c]]$rsid)) %>%
    dplyr::left_join(gwas_hits_pr_rsid[[c]],by=c("rsid")) %>%
    dplyr::mutate(GWAS_HIT = paste0('GWAS_HIT=',GWAS_HIT))
  
  gwas_data_full <- as.data.frame(gwas_hits[[c]] %>%
    dplyr::full_join(dplyr::select(vcf_data, rsid, chrom, 
                                   pos_start, pos_end, 
                                   ref, alt, vartype)) %>%
    dplyr::left_join(gwas_citations[[c]])
  )
  
  fname <- paste0('data/gwas_',datestamp,'.tsv')
  fname_link <- 'gwas.tsv'
  fname_phenotypes <- paste0('data/gwas_phenotypes_',datestamp,'.tsv')
  fname_phenotypes_link <- 'gwas_phenotypes.tsv'
  if(c == 'all'){
    fname <- paste0('data/gwas_all_',datestamp,'.tsv')
    fname_link <- 'gwas_all.tsv'
    fname_phenotypes <- paste0('data/gwas_all_phenotypes_',datestamp,'.tsv')
    fname_phenotypes_link <- 'gwas_all_phenotypes.tsv'
  }
  gwas_phenotypes <- gwas_data_full %>%
    dplyr::mutate(var_id = paste(chrom, pos_start, ref, alt, sep = "_")) %>%
    dplyr::select(var_id, rsid, strongest_snp_risk_allele, 
                  efo_id, p_value, pmid, efo_name, citation, link) %>% 
    dplyr::distinct()
  
  write.table(gwas_phenotypes,file = fname_phenotypes, 
              sep = "\t",col.names = T,row.names = F,quote = F)
  write.table(gwas_data_full, file = fname, 
              sep="\t",col.names = T,row.names = F,quote = F)
  
  system(paste0("ln -F -s ",fname_phenotypes," ",fname_phenotypes_link))
  system(paste0("ln -F -s ",fname," ",fname_link))
  
  
  print_vcf(vcf_data, class = c)
  
}

