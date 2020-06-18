library(tidyverse)
library(reshape2)
library(optparse)

model_path = "../data/mfactorization/credible_set/mapping_sn_spMF_K15_a1800_l1300_Loadings"
variants_metadata_file = "../data/found_variants_in_sumstat/variants_metadata.tsv"
factors_tbl_output =  "../data/found_variants_in_sumstat/eqtl_factors_for_cont_product.tsv"

betas = read.table(paste0(model_path, "_beta.txt"))
pvalues = read.table(paste0(model_path, "_pvalue_BH.txt"))

nrow(betas)
nrow(pvalues)
N_factors = ncol(betas)

colnames(betas) = sapply(strsplit(colnames(betas), split = "V"), "[[", 2)
colnames(pvalues) = sapply(strsplit(colnames(pvalues), split = "X"), "[[", 2)

genes = sapply(strsplit(rownames(betas), split = " "), "[[", 1)
variants = sapply(strsplit(rownames(betas), split = " "), "[[", 2)
eqtl_ids = paste(variants, genes, sep=".")

################# stand beta or max scale ##################
# stand_betas = as.data.frame(t(scale(t(betas), center = F))) # divide by sd
# stand_betas = as.data.frame(t(apply(betas, 1, function(x){abs(x)/max(abs(x))}))) # max scale

# stand_betas$eqtl_id = eqtl_ids
# stand_betas_factors = reshape2::melt(stand_betas, id="eqtl_id", variable.name = "factor", value.name="stand_beta")
############################################################

########################### pip ############################
# cc_file = "../data/connected_components/microarr_cc.tsv"
# ccs = read_tsv(cc_file)
# ccs = ccs %>% mutate(eqtl_id = paste(variant_id, phenotype_id, sep="."))
# 
# pip_df = 
#   ccs %>% 
#   distinct(eqtl_id, cell_type, pip) %>% 
#   group_by(eqtl_id) %>% 
#   summarise(max_pip = max(pip)) 
# 
# stand_betas_factors = expand_grid(pip_df, factor = as.character(seq(1:N_factors)))
# stand_betas_factors = rename(stand_betas_factors, stand_beta = max_pip)
###########################################################

#################### product ##############################
stand_betas = as.data.frame(t(apply(betas, 1, function(x){abs(x)/max(abs(x))})))
stand_betas$eqtl_id = eqtl_ids
stand_betas_factors = reshape2::melt(stand_betas, id="eqtl_id", variable.name = "factor", value.name="stand_beta")

cc_file = "../data/connected_components/microarr_cc.tsv"
ccs = read_tsv(cc_file)
ccs = ccs %>% mutate(eqtl_id = paste(variant_id, phenotype_id, sep="."))

pip_df = 
  ccs %>% 
  distinct(eqtl_id, cell_type, pip) %>% 
  group_by(eqtl_id) %>% 
  summarise(max_pip = max(pip)) 

stand_betas_factors = 
  stand_betas_factors %>% left_join(pip_df, by="eqtl_id") %>% mutate(stand_beta = stand_beta*max_pip)

#########################################################

pvalues$eqtl_id = eqtl_ids
betas$eqtl_id = eqtl_ids

eqtl_factors = reshape2::melt(pvalues, id="eqtl_id", variable.name = "factor", value.name="pvalue")
betas_factors = reshape2::melt(betas, id="eqtl_id", variable.name = "factor", value.name="beta")

eqtl_factors = left_join(eqtl_factors, betas_factors, by=c("eqtl_id","factor")) %>% 
  left_join(stand_betas_factors, by=c("eqtl_id","factor")) 

nrow(eqtl_factors)
eqtl_factors = filter(eqtl_factors, pvalue <= 0.05)
nrow(eqtl_factors)

variants_metadata = read_tsv(variants_metadata_file)
eqtl_factors = left_join(eqtl_factors, variants_metadata, by="eqtl_id")
eqtl_factors %>% distinct(eqtl_id) %>% nrow()

write_tsv(eqtl_factors, factors_tbl_output)

