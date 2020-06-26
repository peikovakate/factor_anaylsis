for chr in {1..22};
do echo ${chr}; 
Rscript factors_to_bed_continuous.R \
  -b  ../data/bim_files/eur_chr${chr}_biall_maf_cm.bim \
  -f ../data/credible_sets/factor_cont_binary.tsv \
  -o ../data/mfactorization/all_variants_from_credible_sets/annot_binary/${chr};
done;