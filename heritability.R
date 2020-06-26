library(tidyverse)
library(ggplot2)

gwas = ""
folder = "../data/mfactorization/all_variants_from_credible_sets/heritability/"
herit = read_tsv(file.path(folder, sprintf("%s.results", gwas)))
nrow(herit)

herit = arrange(herit, Category)

categories = strsplit(herit$Category, split = "L2_0")
categories = unlist(lapply(categories, "[[", 1))

herit = mutate(herit, Category=categories)
significant = herit %>% filter(Enrichment_p <= 0.05) 
# herit = herit %>% filter(Enrichment_p <= 0.05)

ggplot(herit, aes(x=Category, y = Enrichment))+
  ggplot2::geom_col() +
  ggplot2::geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2) + 
  ggplot2::ggtitle(paste("GWAS", toupper(gwas))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ggsave(file.path(folder, paste0(gwas, ".png")))
