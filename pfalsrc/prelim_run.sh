#!/bin/bash

echo -e "Alright let's run some STUFF"
echo -e "Running NFx4026 Specific Curations...\n\n"
python3 NFx4026-specific_curation.py GENE2.count GENE2_NAME.txt NF54_NHP4026_combo.snps.recal.sel.vcf NF54_NHP4026_cursnps.vcf /./data/
echo -e "\n\nNF54x4026 Curation Complete!\n\n"
echo -e "\n\nFixing VCF and Converting GFF to GTF...\n\n"
python3 pfal_genomeandanno_prep.py PlasmoDB-46_Pfalciparum3D7.gff NF54_NHP4026_cursnps.recode.vcf /./data/ 0.9
echo -e "\n\nVCF and GTF Ready!\n\n"
echo -e "\n\nNormalizing count data and generating covariates...\n\n"
mkdir /./data/PEERdata
echo -e "Processing T4..."
python3 normalization_and_covariates.py T4_counts.txt GENE2_probeinfo.csv T4_norm PlasmoDB-46_Pfalciparum3D7.gff.gtf T4_batchcov.csv Gold_Standard.txt NF54_NHP4026_cursnps.recode_annot_final.vcf /./data/ -matrix True
echo -e "Processing T30..."
python3 normalization_and_covariates.py T30_counts.txt GENE2_probeinfo.csv T30_norm PlasmoDB-46_Pfalciparum3D7.gff.gtf T30_batchcov.csv Gold_Standard.txt NF54_NHP4026_cursnps.recode_annot_final.vcf /./data/ -matrix True
echo -e "Processing T44..."
python3 normalization_and_covariates.py T44_counts.txt GENE2_probeinfo.csv T44_norm PlasmoDB-46_Pfalciparum3D7.gff.gtf T44_batchcov.csv Gold_Standard.txt NF54_NHP4026_cursnps.recode_annot_final.vcf /./data/ -matrix True

echo -e "\n\nNormalization and covariate generation complete!\n\n"

echo -e "Running PEER..."
for i in 1 2 3 4 5 10 15 20 25 30
do
	Rscript run_PEER.R /./data/T4_norm.bed.gz T4_norm_$i $i --max_iter 1000 -o /data/PEERdata -c /data/T4_batchcov_crc.txt
	Rscript run_PEER.R /./data/T30_norm.bed.gz T30_norm_$i $i --max_iter 1000 -o /data/PEERdata -c /data/T30_batchcov_crc.txt
	Rscript run_PEER.R /./data/T44_norm.bed.gz T44_norm_$i $i --max_iter 1000 -o /data/PEERdata -c /data/T44_batchcov_crc.txt
done
echo "PEER finished!"