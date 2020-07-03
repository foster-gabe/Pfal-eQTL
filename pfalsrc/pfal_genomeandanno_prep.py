#!/usr/bin/env python
# coding: utf-8


################################################################
#
# Title: eQTL File Preparation
# Author: Gabe Foster
# Date: 6/5/2020
# Purpose: This script prepares the various files necessary
# for an eQTL run.
################################################################

import pandas as pd
import os
import argparse
import re



def gff_to_gtf(gff_file):

    '''
	This is a quick function to convert the PlasmoDB gff file for
	P. falciparum in to a gtf format with gene level information
	only
    '''

    gff_file = pd.read_csv(gff_file,
                           sep = '\t',
                           comment = '#',
                           header = None)

    gtf_data = []
    gff_file = gff_file[gff_file.iloc[:,2] == 'gene']
    gff_file.iloc[:,8] = gff_file.iloc[:,8].str.replace('^ID=', '', regex = True)
    gff_file.iloc[:,8] = gff_file.iloc[:,8].str.replace(';.*', '', regex = True)
    for iter, row in gff_file.iterrows():
        row[8] = f'gene_id "{row[8]}";'
        gtf_data.append(row)
    gtf_out = pd.DataFrame(gtf_data)
    return gtf_out


def prepare_vcf(vcf, path, min_snp_coverage):
    
    '''
    This piece prepares vcf files for cis and trans analysis- the trans vcf file
    needs the markers to be pruned for LD. Additionally, we need some form of variant ID in here;
    this is just going to be chr~geneloc. It does work.
    This piece requires the 'plink_pruning_prep.sh script.
    '''
	
	# Calls a little bash script that juggles things to allow PLINK to do the heavy lifting-
	# removes bad markers, removes markers with sig LD
	
    missing_snp_pct = 1 - min_snp_coverage
	
    vcfhandle = vcf.replace('.vcf', '')
	
    os.system(f'plink --indep 50 5 2 --id-delim "~" --vcf {vcf} --allow-extra-chr --const-fid A --out {path}tmp1 --geno {missing_snp_pct}')
    os.system(f'plink --vcf {vcf} --extract {path}tmp1.prune.in --make-bed --out {path}tmp2 --const-fid A --allow-extra-chr') 
    os.system(f'plink --bfile {path}tmp2 --recode vcf-iid --out {vcfhandle}_final --allow-extra-chr')
	
    os.system(f'./plink_pruning_prep.sh {vcf} {missing_snp_pct}')
    os.system(f'cp *_final.vcf {path}')
	
    os.system('rm /./data/tmp2* /./data/tmp1*')

		
	# These lines will gzip and index the vcf; we are curating the vcfs later
	# for our use, so these are not in play now
	
	# Easy, now we just zip it up and index it
	
    # os.system(f'bgzip {path}{vcfhandle}_final.vcf')
    # os.system(f'tabix -p vcf {path}{vcfhandle}_final.vcf.gz')
    
    


def to_gct(file, expdata):
	mrow = 0
	mcol = 0
	row = len(expdata.index)
	col = len(expdata.columns)
	outfile = open(file, 'w')
	outfile.write(f'1.0\n')
	outfile.write(f'{row}\t{col}\n')
	expdata.to_csv(file,
                       mode = 'a',
                       header = True,
                       sep = '\t',
                       encoding='utf-8',
                       index = True,
					   index_label = False)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Reformats the PlasmoDB GFF to a GTF and Indexes the vcf file')
    parser.add_argument('gff_file', help = 'The PlasmoDB .GFF file')
    parser.add_argument('vcf_file', help = 'the *.vcf file from this genetic cross')
    parser.add_argument('data_path', help = 'path to data folder')
    parser.add_argument('min_snp_coverage', help = 'minimum percentage of samples with a snp call for each marker, pct', default = 0.90, type = float)
    args = parser.parse_args()
    
    if (args.gff_file):
        gtf = gff_to_gtf(f'{args.data_path}{args.gff_file}')
        gtf_handle = re.sub('r[\.[^.]+$','',args.gff_file)
        gtf_name = (f'{args.data_path}{gtf_handle}.gtf')
        gtf.to_csv(gtf_name, sep = '\t', index = False, header = False)
    else:
        print('No GFF file provided')
    
    if (args.vcf_file):
        prepare_vcf(f'{args.data_path}{args.vcf_file}', args.data_path, args.min_snp_coverage)
    else:
        print('No VCF provided')
    
  

