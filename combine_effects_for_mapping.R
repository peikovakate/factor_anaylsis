library(tidyverse)
library(reshape2)

variants_file = "../data/sumstat_cs_variants/tissues_combined.tsv"
metadata_file = "../data/sumstat_cs_variants/variants_metadata.tsv"
tissues_file = "../data/sumstat_cs_variants/tissues.txt"
output_dir = "../data/sumstat_cs_variants/"

# prepare data for mapping (all) eqtls to factorization matrix
found_variants = read_tsv(variants_file)

# if cell name contains directory as well
if(grepl("/", found_variants$cell_type[1])){
  cell_types = strsplit(found_variants$cell_type, split = "/")
  cell_types = unlist(lapply(cell_types, "[[", 2)) 
  found_variants = mutate(found_variants, cell_type = cell_types)
}

# if columns got renamed because of query_sumstat_for_variants.R script
# rename them back
if("variant_id" %in% colnames(found_variants)){
  found_variants = found_variants %>% 
    rename(variant = variant_id, molecular_trait_id = phenotype_id)  
}

found_variants = found_variants %>% 
  mutate(eqtl_id = paste(variant, molecular_trait_id, sep="."))

# remove identical records
found_variants = found_variants %>% distinct(.keep_all = T)

metadata = found_variants %>% 
  select(eqtl_id, molecular_trait_id, variant, chromosome, position) %>% 
  distinct(.keep_all = T) 

write_tsv(metadata, metadata_file)

eqtls <- found_variants %>%
  group_by(eqtl_id) %>%
  mutate(count = length(unique(cell_type))) %>% 
  ungroup()

eqtls %>%
  group_by(eqtl_id,cell_type) %>%
  filter (pvalue == min(pvalue)) %>% 
  ungroup() %>% 
  nrow()

eqtls %>% 
  distinct(eqtl_id, cell_type) %>% 
  nrow()

eqtls_present_in_all_tissues = filter(eqtls, count == 17)

eqtls_present_in_all_tissues %>% 
  distinct(eqtl_id, cell_type) %>% 
  nrow()

pvalues = dcast(eqtls_present_in_all_tissues, eqtl_id ~ cell_type, value.var = "pvalue") %>% dplyr::arrange(eqtl_id)
betas = dcast(eqtls_present_in_all_tissues, eqtl_id ~ cell_type, value.var = "beta") %>% dplyr::arrange(eqtl_id)
ses = dcast(eqtls_present_in_all_tissues, eqtl_id ~ cell_type, value.var = "se") %>% dplyr::arrange(eqtl_id)

write_lines(setdiff(colnames(betas), c("eqtl_id")), tissues_file)

# write_tsv(pvalues[, 2:ncol(pvalues)], "../data/found_variants_in_sumstat/pvalues.txt")
# write_tsv(betas[, 2:ncol(betas)], "../data/found_variants_in_sumstat/betas.txt")
# write_tsv(ses[, 2:ncol(ses)], "../data/found_variants_in_sumstat/ses.txt")

# columns = c("eqtl_id", "molecular_trait_id", "chromosome", "position", "variant", "cc_id", "count")
columns = c("eqtl_id", "molecular_trait_id", "variant")
eqtl_info = distinct(eqtls_present_in_all_tissues[, columns], .keep_all = T)
pvalues = left_join(pvalues, eqtl_info, by="eqtl_id")
betas = left_join(betas, eqtl_info, by="eqtl_id")
ses = left_join(ses, eqtl_info, by="eqtl_id")

# write_tsv(pvalues, "../data/found_variants_in_sumstat/eqtl_pval_present_in_all_tissues.tsv")
# write_tsv(betas, "../data/found_variants_in_sumstat/eqtl_beta_present_in_all_tissues.tsv")
# write_tsv(ses, "../data/found_variants_in_sumstat/eqtl_se_present_in_all_tissues.tsv")
# write_tsv(ses, "../data/found_variants_in_sumstat/eqtl_se_present_in_all_tissues.tsv")

betas %>% 
  rename(SNP = variant, Gene = molecular_trait_id) %>% 
  select(-eqtl_id) %>% 
  select(SNP, Gene, everything()) %>%  
  write_tsv(file.path(output_dir, "slope.txt"))

ses %>% 
  rename(SNP = variant, Gene = molecular_trait_id) %>% 
  select(-eqtl_id) %>% 
  select(SNP, Gene, everything()) %>% 
  write_tsv(file.path(output_dir, "se.txt"))



