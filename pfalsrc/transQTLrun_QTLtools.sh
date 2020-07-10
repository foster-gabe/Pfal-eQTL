#!/bin/sh

#Let's run some QTLLLLLLLLLLLLLL

mkdir transeQTL

# For each timepoint run the entire analysis:

for k in T4 T44
do
	# For each # inferred covariates being tested:
	for j in 3
	do
	
		QTLtools trans --vcf $k\_norm.vcf.gz --bed $k\_norm.bed.gz --cov PEERdata/$k\_norm_$j.PEER_covariates.txt --out transeQTL/$k\_transperm_$j --normal --sample 1000 --window 1000000
		gunzip transeQTL/$k\_transperm_$j.best.txt.gz
		cat transeQTL/$k\_transperm_$j.best.txt | cut -d" " -f1-3 > transeQTL/$k\_transperm_$j.best_cut.txt
		gzip transeQTL/$k\_transperm_$j.best_cut.txt
		QTLtools trans --vcf $k\_norm.vcf.gz --bed $k\_norm.bed.gz --cov PEERdata/$k\_norm_$j.PEER_covariates.txt --adjust transeQTL/$k\_transperm_$j.best_cut.txt.gz --out transeQTL/$k\_transadjust_$j --normal --threshold .5 --window 1000000
		
	# For each significance level:
	
		for i in 0.05 0.10 0.25
		do
			
			Rscript runFDR_atrans.R transeQTL/$k\_transadjust_$j.best.txt.gz transeQTL/$k\_transadjust_$j.hits.txt.gz $i transeQTL/$k\_trans_FDR_$j\_$i.txt
	
		done
	done
done

for k in T30
do
	# For each # inferred covariates being tested:
	for j in 4
	do
	
		QTLtools trans --vcf $k\_norm.vcf.gz --bed $k\_norm.bed.gz --cov PEERdata/$k\_norm_$j.PEER_covariates.txt --out transeQTL/$k\_transperm_$j --normal --sample 1000 --window 1000000
		gunzip transeQTL/$k\_transperm_$j.best.txt.gz
		cat transeQTL/$k\_transperm_$j.best.txt | cut -d" " -f1-3 > transeQTL/$k\_transperm_$j.best_cut.txt
		gzip transeQTL/$k\_transperm_$j.best_cut.txt
		QTLtools trans --vcf $k\_norm.vcf.gz --bed $k\_norm.bed.gz --cov PEERdata/$k\_norm_$j.PEER_covariates.txt --adjust transeQTL/$k\_transperm_$j.best_cut.txt.gz --out transeQTL/$k\_transadjust_$j --normal --threshold .5 --window 1000000
		
	# For each significance level:
	
		for i in 0.05 0.10 0.25
		do
			
			Rscript runFDR_atrans.R transeQTL/$k\_transadjust_$j.best.txt.gz transeQTL/$k\_transadjust_$j.hits.txt.gz $i transeQTL/$k\_trans_FDR_$j\_$i.txt
	
		done
	done
done