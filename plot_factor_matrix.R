library(pheatmap)

matrix_file = "../data/rnaseq/unique_genes/output_rnaseq/sn_spMF_K35_a1600_l11000_Run4.RData"
load(matrix_file)
pheatmap(FactorM)

