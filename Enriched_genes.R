library(readr)
library(readxl)
library(easyPubMed)
library(multidplyr)
library('org.Hs.eg.db')
library('org.Rn.eg.db')
library('org.Mm.eg.db')
library(tidyverse)
library(RISmed)
library(future.apply)
library(writexl)

# get the mapping between human ENTREZ ids and gene symbols
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
all((as.list(x[mapped_genes]) %>% map_int(length)) == 1)
human_g <- enframe(as.list(x[mapped_genes])) %>% mutate(value=unlist(value))

# get the mapping between rat ENTREZ ids and gene symbols
x <- org.Rn.egSYMBOL
mapped_genes <- mappedkeys(x)
all((as.list(x[mapped_genes]) %>% map_int(length)) == 1)
rat_g <- enframe(as.list(x[mapped_genes])) %>% mutate(value=unlist(value))

# get the mapping between mouse ENTREZ ids and gene symbols
x <- org.Mm.egSYMBOL
mapped_genes <- mappedkeys(x)
all((as.list(x[mapped_genes]) %>% map_int(length)) == 1)
mouse_g <- enframe(as.list(x[mapped_genes])) %>% mutate(value=unlist(value))

# load the table with the mapping between SARS-CoV-2 ENTREZ ids and SARS-CoV-2 gene symbols
sars_genes <- read_delim("sars_genes.txt", 
                         "\t", escape_double = FALSE, col_types = cols(GeneID = col_character()), 
                         trim_ws = TRUE) %>% dplyr::select(name=GeneID,value=Symbol)

mapping_tab <- bind_rows(human_g=human_g, rat_g=rat_g, mouse_g=mouse_g, sars_genes=sars_genes,.id = "species") %>% 
              dplyr::select(GeneID=name,symbol=value,everything())

#load the table PMC-ids.csv table ( https://ftp.ncbi.nlm.nih.gov/pub/pmc/PMC-ids.csv.gz )
PMC_ids <- read_csv("PMC-ids.csv", col_types = cols(PMID = col_character()))

#load the gene2pubmed table ( https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz )
gene2pubmed <- read_delim("gene2pubmed","\t", escape_double = FALSE, 
                          col_types = cols(`#tax_id` = col_character(), 
                                           GeneID = col_character(), PubMed_ID = col_character()), 
                          trim_ws = TRUE)

#load the table listing COVID19 papers used in this manuscript
docs_filtered_new <- read_excel("table-27k.xlsx")

# Associate the publication year to each entry of the gene2pubmed table 
# (needed as a prefilter to reduce the number of manuscripts when querying NCBI)
anno <- gene2pubmed %>% left_join(dplyr::select(PMC_ids,"PMID","DOI","Year"), by = c("PubMed_ID" = "PMID"))

#keep on associations involving human/rat/mouse/Sars genes and manuscripts published in 2020
anno_20 <- anno %>% filter(GeneID %in% mapping_tab$GeneID, Year %in% c(2020))

#Annotate manuscripts 
anno_group <- anno_20 %>% group_by(PubMed_ID,Year) %>% summarize(n_genes = n()) %>% ungroup()
anno_group$is_covid <- anno_group$PubMed_ID %in% docs_filtered_new$accn

# function to retrieve more detailed information from Pubmed 
get_info <- function(pmid){
  my_query <- paste0(pmid,'[PMID]')
  my_entrez_id <- get_pubmed_ids(my_query)
  my_abstracts_txt <- fetch_pubmed_data(my_entrez_id)
  df <- article_to_df(my_abstracts_txt,getAuthors = F, getKeywords=T)
  res <- EUtilsSummary(my_query)
  df$year <- YearPubmed(EUtilsGet(res))
  df$month <- MonthPubmed(EUtilsGet(res))
  df$day <- DayPubmed(EUtilsGet(res))
  df
}

#function to manually set the date of manuscripts for which automatic annoataion ha failed
manSetDate <- function(df,PubMed_ID,year,month, day ){
  df[which(df$PubMed_ID == PubMed_ID),]$year <- year
  df[which(df$PubMed_ID == PubMed_ID),]$month <- month
  df[which(df$PubMed_ID == PubMed_ID),]$day <- day
  df
}


# fro each paper in anno_group retrieve information from pubmed
dois <- list()
for (i in 1:length(anno_group$PubMed_ID)){
  print(i)
  dois[[ anno_group$PubMed_ID[i] ]] <- tryCatch(get_info(anno_group$PubMed_ID[i]), error=function(err) NA)
}
#build a table with needed annotation info for each manuscript 
dois_df <- dois %>% do.call(rbind,.)
dois_df <- dois_df[,1:10] %>% as_tibble()
anno_group <- anno_group %>% left_join(dois_df, by=c("PubMed_ID" = "pmid"))

#filter annotated paper to match the same time frame of those reported in the table listing COVID19 papers used in this manuscript
anno_group <- anno_group %>% dplyr::filter(year == 2020 & ((month == 6 & day <= 17) | month < 6 ))

#get the association gene/manuascript for all papers annotated in anno_group
anno_20_filt <- anno_20 %>% filter(PubMed_ID %in% anno_group$PubMed_ID)

#count how many times each gene is cited in the whole dataset
anno_group_gene <- anno_20_filt %>% left_join(anno_group,by=("PubMed_ID")) %>% group_by(GeneID) %>% dplyr::summarise(n=n())

#count how many times each gene is cited in the covid dataset (the list of papers)
anno_group_gene_covid <- anno_20_filt %>% left_join(anno_group,by=("PubMed_ID")) %>% filter(is_covid) %>% 
                         group_by(GeneID) %>% dplyr::summarise(covid_papers=n())

#count how many times each gene is cited in the non-covid dataset ( 2020 papers non in the list)
anno_group_gene_not_covid <- anno_20_filt %>% left_join(anno_group,by=("PubMed_ID")) %>% filter(!is_covid) %>% 
                        group_by(GeneID) %>% dplyr::summarise(not_covid_papers=n())

# merge information for each geneID 
anno_group_gene_cov <- anno_group_gene %>% left_join(anno_group_gene_covid) %>% left_join(anno_group_gene_not_covid) %>% 
                        mutate(covid_papers = if_else(is.na(covid_papers),0L,covid_papers)) %>% 
                        mutate(not_covid_papers = if_else(is.na(not_covid_papers),0L,not_covid_papers))

# Summarize information at the gene symbol level
anno_group_gene_cov <- anno_group_gene_cov %>% 
                        left_join(mapping_tab, by="GeneID") %>%
                        mutate(symbol=toupper(symbol)) %>% 
                        group_by(symbol) %>% 
                        summarize(n = sum(n),
                                  covid_papers = sum(covid_papers), 
                                  not_covid_papers = sum(not_covid_papers)) %>% 
                        ungroup()
  
# get total cont of covid and not-covid manuscripts
univ <- table(anno_group$is_covid)

plan(multisession, workers = 100)

#Compute hypergeometric test for each gene (test if it is enriched in one of the two categories of manuscripts)
anno_group_gene_cov$pvals <- future_apply(anno_group_gene_cov, 1, function(x) {
                                  tbl <- matrix(as.numeric(c(x[4:3],univ)), ncol=2, byrow=T)
                                  fisher.test(tbl, alternative="two.sided")$p.value
                                  })

#adjust p-values for multiple testing
anno_group_gene_cov$pvals_fdr <- p.adjust(anno_group_gene_cov$pvals, method="fdr")
anno_group_gene_cov <- anno_group_gene_cov %>% arrange(pvals_fdr)

#write reults to .xlsx file
write_xlsx(anno_group_gene_cov, "Tab_2020_cov_fisher_h_m_r.xlsx")

#save image
save.image("covid_lit_20.Rdata")
