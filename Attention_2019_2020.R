library(readr)
library(readxl)
library(tidyverse)
library(multidplyr)
library('org.Hs.eg.db')
library('org.Rn.eg.db')
library('org.Mm.eg.db')
library(writexl)
library(furrr)
library(purrr)

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

# Associate the publication year to each entry of the gene2pubmed table 
anno <- gene2pubmed %>% left_join(dplyr::select(PMC_ids,"PMID","DOI","Year"), by = c("PubMed_ID" = "PMID"))

#keep on associations involving human/rat/mouse/Sars genes and manuscripts published in 2019 and 2020
anno_19_20 <- anno %>% filter(Year %in% c(2019,2020), GeneID %in% mapping_tab$GeneID)

anno_group <- anno_19_20 %>% group_by(PubMed_ID,Year) %>% dplyr::summarize(n_genes = n()) %>% ungroup()

#the data table is a gene-manuscript association table containing for each gene-manuscript association 
# the geneID, gene symbol and PubMed_ID, Year
# this function ranks genes in data table based on the number of times they have been cited in each year
# then computes (for each gene that is cited more than 5 times in 2020) the rank difference between 2019 and 2020 

rank_jump <- function(data) {
  #Annotate manuscripts from the gene/manuscript 
  anno_group <- data %>% group_by(PubMed_ID,Year) %>% dplyr::summarize(n_genes = n()) %>% ungroup()
  anno_group$year_2020 <- anno_group$Year == 2020

  #count how many times each gene is cited in the whole dataset
  anno_group_gene <- data %>% left_join(anno_group,by=("PubMed_ID")) %>% group_by(GeneID) %>% dplyr::summarise(n=n())
  
  #count how many times each gene is cited in the 2020
  anno_group_gene_2020 <- data %>% left_join(anno_group,by=("PubMed_ID")) %>% 
                                   filter(year_2020) %>% 
                                   group_by(GeneID) %>% 
                                   dplyr::summarise(papers_2020=n())
  
  #count how many times each gene is cited in the 2019
  anno_group_gene_2019 <- data %>% left_join(anno_group,by=("PubMed_ID")) %>% 
                                   filter(!year_2020) %>% 
                                   group_by(GeneID) %>% 
                                   dplyr::summarise(papers_2019=n())
  
  ## merge information for each geneID 
  anno_group_gene_cov <- anno_group_gene %>% left_join(anno_group_gene_2020) %>% left_join(anno_group_gene_2019) %>% 
                                            mutate(papers_2020 = if_else(is.na(papers_2020),0L,papers_2020)) %>% 
                                            mutate(papers_2019 = if_else(is.na(papers_2019),0L,papers_2019))
  
  # Summarize information at the gene symbol level
  anno_group_gene_cov <- anno_group_gene_cov %>% 
    left_join(mapping_tab, by="GeneID") %>%
    mutate(symbol=toupper(symbol)) %>% 
    group_by(symbol) %>% 
    summarize(n = sum(n),
              papers_2020 = sum(papers_2020), 
              papers_2019 = sum(papers_2019)) %>% 
    ungroup()
  
  #rank gene symbols based on the number of citations in each year and compute the difference
  ranked_genes <- anno_group_gene_cov %>% 
                dplyr::filter(papers_2020 >= 5) %>% 
                mutate( rank_19 = dplyr::min_rank(dplyr::desc(papers_2019)),
                        rank_20 = dplyr::min_rank(dplyr::desc(papers_2020)),
                        rank_jump = -(rank_20 - rank_19))
  
  ranked_genes[,c("symbol","rank_jump")]
}


# produce a random realization of the gene association table
# by sampling an equal number of papers from each year and 
# permuting genes on the association tables in each realization

bootstrap <- function(x, data){

  min_sam <- min(table(data$Year))
  
  #sample an equal number of papers from each year
  samp_2019 <- data %>% filter(Year == 2019) %>% sample_n(min_sam)
  samp_2020 <- data %>% filter(Year == 2020) %>% sample_n(min_sam)
  #permute genes in the association table
  perm <- transform( bind_rows(samp_2019,samp_2020), GeneID = sample(GeneID,) )
  rank_jump(perm)
}

# Set a "plan" for how the code should run.
plan(multisession, workers = 100)

#perform 1000 iteration for the bootstrap procedure
boot_out <- future_map(1:1000, bootstrap, anno_19_20, .progress = T)

#join the results from each random realization
boot_out_df <- boot_out %>% purrr::reduce(full_join, by = "symbol")
colnames(boot_out_df)[2:ncol(boot_out_df)] <- paste0("sim",1:1000)

#compute the observed differences in ranks between citations from 2019 and from 2020
observed <- rank_jump(anno_19_20) %>% transmute(symbol, observed_jump=rank_jump)

boot_out_df <- observed %>% full_join(boot_out_df, by = "symbol")

#compute empiric pvalues
boot_out_df$p <- unlist(apply(boot_out_df,1, function(x){
                                                  sum(x[3:length(x)] >= x[2],na.rm = T)/1000
                                              }
                              )
                        ) 

#filter the complete datatset for the genes present in the observed dataset
boot_out_df_sym <- boot_out_df %>% dplyr::select(symbol, observed_jump, p, everything()) %>% 
                                   dplyr::filter(!is.na(observed_jump))

# make a selection of top genes 
sig <- boot_out_df_sym %>% filter(p<0.10, observed_jump > 0) %>% 
                           arrange(desc(observed_jump)) %>% 
                           dplyr::select(symbol, observed_jump, p)


saveRDS(boot_out_df_sym, "boot_out.RDS")

writexl::write_xlsx(boot_out_df, "bootstrap.xlsx")
writexl::write_xlsx(sig, "sig_bootstrap.xlsx")

