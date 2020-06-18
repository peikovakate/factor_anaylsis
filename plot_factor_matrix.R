library(pheatmap)

matrix_file = "../data/mfactorization/output/sn_spMF_K15_a1800_l1300/sn_spMF_K15_a1800_l1300_Run28.RData"
load(matrix_file)
pheatmap(FactorM)

