#!/bin/sh

#Let's run some QTLLLLLLLLLLLLLL

mkdir ciseQTLcovariate
j=4
for j in 1 2 3 4 5 10 15 20 25 30 nobatch nocrc
do
	for i in $(seq 20)
	do
	
	fastQTL -V T4_norm.vcf.gz -B old_T4_norm.bed.gz --permute 1000  -C PEERdata/T4_norm_$j.PEER_covariates.txt --out ciseQTLcovariate/T4_cis_$j\_$i.txt --chunk $i 20  

	done

	for i in $(seq 20)
	do

	fastQTL -V T30_norm.vcf.gz -B old_T30_norm.bed.gz --permute 1000  -C PEERdata/T30_norm_$j.PEER_covariates.txt --out ciseQTLcovariate/T30_cis_$j\_$i.txt --chunk $i 20  

	done

	for i in $(seq 20)
	do

	fastQTL -V T44_norm.vcf.gz -B old_T44_norm.bed.gz --permute 1000  -C PEERdata/T44_norm_$j.PEER_covariates.txt --out ciseQTLcovariate/T44_cis_$j\_$i.txt --chunk $i 20  

	done

done

