#!/bin/sh

#Let's run some QTLLLLLLLLLLLLLL

mkdir ciseQTL

# For each timepoint run the entire analysis:

for k in T4 T30 T44
do
	# For each # inferred covariates being tested:
	for j in 0 1 2 3 4 5 10 15 20 25 30 nobatch nocrc
	do
	
		QTLtools cis --vcf $k\_norm.vcf.gz --bed $k\_norm.bed.gz --cov PEERdata/$k\_norm_$j.PEER_covariates.txt --out ciseQTL/$k\_cisnom_$j.txt --nominal 1 --window 1000000
		QTLtools cis --vcf $k\_norm.vcf.gz --bed $k\_norm.bed.gz --cov PEERdata/$k\_norm_$j.PEER_covariates.txt --out ciseQTL/$k\_cisperm_$j.txt --permute 1000 10000 --window 1000000
	
	# For each significance level:
	
		for i in 0.05 0.10 0.25
		do
			Rscript runFDR_cis.R ciseQTL/$k\_cisperm_$j.txt $i ciseQTL/$k\_FDRcuts_$i\_$j.txt
			Rscript NomFilter.R ciseQTL/$k\_cisnom_$j.txt ciseQTL/$k\_FDRcuts_$i\_$j.txt.thresholds.txt ciseQTL/$k\_cisnom_sig_$j.txt
	
		done
	done
done

