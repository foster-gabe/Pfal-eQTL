#!/usr/bin/env python
# coding: utf-8


######################################################################
#
# Name: Normalization and Covariate Generation
# Author: Gabe Foster
# Date: 6/7/2020
# Purpose: This script performs TPM normalization on count data in
# gct format, handles staging by transcription and generation of
# cyclical regression covariates and prepares files for PEER
#
######################################################################

import pandas as pd
import pfal_genomeandanno_prep as pfal
import os
import argparse
import numpy as np

def sample_and_geneqc(expdata):
    '''
    Here, we'll remove poor samples and poor genes, and QC low read #s
    Low reads are any gene/sample reading < 5 counts; those are zeroed out
    Poor samples are defined as those with < 3000 genes with reads
    Poor genes are those that appear in < 20% of samples
    '''
    
    # Curate Samples
    
    expdata = expdata.mask(expdata < 5, 0)
    
    genecounts = pd.Series(data = np.count_nonzero(expdata, axis = 1),
                             index = expdata.index)
    samplecounts = pd.Series(np.count_nonzero(expdata,axis = 0),
                          index = expdata.columns)
    
    goodgenes = genecounts[genecounts/samplecounts.size > 0.2]
    
    goodsamples = samplecounts[samplecounts > 3000]
    
    allcur = expdata.loc[goodgenes.index, goodsamples.index]
    
    return allcur
	
def tpm_norm(expdata, probefile, fraglength = 250):
    '''
    TPM is (reads / effective gene length (in kb)) / (sample reads / 1e6)
    The elaborate nonsense below just calculates this, it's a bit of a
    mess but hey that's programming
    '''
    
    # Read in probe metadata, subtract frag length,
    # and convert to dict
    
    
    probeinfo = pd.read_csv(probefile, usecols = [1,6])
    probeinfo.index = probeinfo['Geneid']
    probeinfo.drop(columns = 'Geneid', inplace = True)
    probeinfo = probeinfo[probeinfo.index.isin(gctdata.index)]
    probeinfo['Length'] = probeinfo['Length'].apply(lambda x: x-fraglength)
        
    # Build output frame from input frame (cheating), calc length and
    # build lookup from samplename to total counts
    
    tpm = expdata.copy()
    curated_tpm = expdata.copy()
    
    # Iterate over rows- if non sample rows, just copy data, if sample
    # rows, calculate TPM
    
    # Note: it really doesn't like the way I did this and will throw
    # warnings- it does work, it's fine, don't worry about it, if you
    # are smarter than I am go ahead and rewrite it
    
    for (col, data) in expdata.iteritems():
        
        numreads = data.sum()
        
        tempframe = pd.concat([probeinfo, data], axis = 1)

        tempframe['rpkb'] = tempframe.iloc[:,1].divide(tempframe['Length'].divide(1000))
        tempframe['tpm'] = (tempframe['rpkb'] / data.sum()) * 1e6
        tempframe['tpm'] = tempframe['tpm'].clip(lower = 0)
        curated_tpm.loc[:,col] = tempframe['tpm']
    	
    # Removing genes whose expression is ~0 in > 80% of samples
    goodgenes=[]
    for index, row in curated_tpm.iterrows():
	
        if (np.count_nonzero(row > 0.1) / len(row) > 0.2):
            goodgenes.append(index)
    	
    curated_tpm = curated_tpm[curated_tpm.index.isin(goodgenes)]
    	
		
    return curated_tpm
    

def create_bed(expdata, gtf_file):
    '''
    Time to create our bed file for QTL analysis
    Note- this is a specific BED format for QTLtools
    '''
    
    # Read file
    gtfdata = pd.read_csv(gtf_file, sep = '\t', index_col = False)
    
    # Fix gene names; current format is '"gene_id PF...;"
    gtfdata.columns = ['chr', 'source', 'type', 'start', 'end','period', 'strand', 'period2', 'geneid']
    gtfdata['geneid'] = gtfdata['geneid'].str.replace('^gene_id', '', regex = True)
    gtfdata['geneid'] = gtfdata['geneid'].str.replace('"', '', regex = True)
    gtfdata['geneid'] = gtfdata['geneid'].str.replace(';', '', regex = True)
    gtfdata['geneid'] = gtfdata['geneid'].str.replace(' ', '', regex = True)
    
    # Pull columns required by QTLtools BED format, sort
    
    beddata = gtfdata.loc[:,['chr','start','end','geneid','geneid','strand']]
    beddata.columns.values[3] = 'exonid'
    beddata.sort_values(by=['chr', 'start'], inplace = True)
    
    # Merge with exp data on geneid, drop redundant cols, add comment char to first row
    
    full_bed = beddata.merge(expdata, left_on = 'exonid', right_on = 'Geneid', how = 'inner')

    full_bed.columns.values[0] = '#chr'   
    
    return full_bed
	
def fqtl_bed(expdata, gtf_file):
    '''
    So QTLtools won't run cis, and fastQTL won't run trans, so we'll use both.
    This function creates a fastQTL compatible BED file.
    '''
    # Read file
    gtfdata = pd.read_csv(gtf_file, sep = '\t', index_col = False)
    
    # Fix gene names; current format is '"gene_id PF...;"
    gtfdata.columns = ['chr', 'source', 'type', 'start', 'end','period', 'strand', 'period2', 'ID']
    gtfdata['ID'] = gtfdata['ID'].str.replace('^gene_id', '', regex = True)
    gtfdata['ID'] = gtfdata['ID'].str.replace('"', '', regex = True)
    gtfdata['ID'] = gtfdata['ID'].str.replace(';', '', regex = True)
    gtfdata['ID'] = gtfdata['ID'].str.replace(' ', '', regex = True)
    beddata = gtfdata.loc[:,['chr','start','end','ID']]
    beddata.sort_values(by=['chr', 'start'], inplace = True)
    full_bed = beddata.merge(expdata, left_on = 'ID', right_on = 'Geneid', how = 'inner')

    full_bed.columns.values[0] = '#chr'   
    
    return full_bed

def index_bed(bed_file):
    '''
    This function just runs the simple tabix call to index the BED file
    '''
    bed_prefix = bed_file.split('\.')[0]
    
    os.system(f'bgzip {bed_file} && tabix -p bed {bed_prefix}.gz')
	
def curate_vcf(samples, mainvcf, alias_out, data_path, matrix):
    '''
    This function creates a vcf file that only contains the samples in the specific set;
    this is necessary for fastQTL trans analysis.
    '''
	
    os.system(f'bcftools view -Oz -s {samples} -o {data_path}{alias_out}.vcf.gz {data_path}{mainvcf}')
    os.system(f'tabix -p vcf {data_path}{alias_out}.vcf.gz')
    
	# Some programs (*cough* MatrixEQTL *cough*) don't take a vcf, and require a matrix instead;
	# this piece creates that by using the vcftools -012 function and combining the matrix, samples
	# and marker files in to one named matrix.
	
    if matrix == True:   	
    	os.system(f'vcftools --vcf {data_path}{alias_out}.vcf --012 --out {data_path}{alias_out}')
    	matrix = pd.read_csv(f'{data_path}{alias_out}.012', sep = '\t', header = None, index_col = 0)
    	sample_names = pd.read_csv(f'{data_path}{alias_out}.012.indv', sep = '\t', header = None, index_col = False)
    	sample_names.columns = ['samples']
    	marker_names = pd.read_csv(f'{data_path}{alias_out}.012.pos', sep = '\t', header = None, index_col = False)
    	marker_names.columns = ['chr', 'pos']
    	marker_names['pos'] = marker_names['pos'].astype(str)
    	marker_names['ID'] = marker_names['chr'].str.cat(marker_names['pos'], sep = '~')
    	matrix.index = sample_names['samples']
    	matrix.columns = marker_names['ID']
    	matrix.to_csv(f'{data_path}{alias_out}.matrix', sep = '\t')
	
 


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Performs TPM normalization and prepares covariates for PEER')
    parser.add_argument('exp_file', help = 'QCed count data as tab delimited text file, NxG')
    parser.add_argument('probe_info', help = 'probe information file from annotation preparation')
    parser.add_argument('alias_out', help = 'file alias for normalized output')
    parser.add_argument('-fraglength', help = 'mean fragment read length, default is 250bp', default = 250)
    parser.add_argument('gtf_file', help = 'GTF file for BED conversion')
    parser.add_argument('cov_file', help = 'Covariate File (csv)')
    parser.add_argument('ref_course', help = 'Reference Time Course (currently Gold_Standard.txt)')
    parser.add_argument('vcffile', help = 'VCF file for this cross')
    parser.add_argument('data_path', help = 'Path to data files')
    parser.add_argument('-matrix', help = 'add genotype matrix', default = False)
    args = parser.parse_args()
	
	    
    # Read in expression data
    
    gctdata = pd.read_csv(f'{args.data_path}{args.exp_file}', 
                          sep = '\t',
    					  index_col = 0)
	
    qcexp = sample_and_geneqc(gctdata)
    
    normcount = tpm_norm(qcexp, f'{args.data_path}{args.probe_info}')
    normcount.to_csv(f'{args.data_path}{args.alias_out}.txt',
    index = True,
    header = True,
    sep = '\t')
	
    samples = list(normcount.columns)
    samples = ','.join(samples)
    curate_vcf(samples, args.vcffile, args.alias_out, args.data_path, args.matrix)
    
    out_name = f'{args.alias_out}.bed'
    beddata = create_bed(normcount, f'{args.data_path}{args.gtf_file}')
    beddata.to_csv(f'{args.data_path}{out_name}', sep = '\t', index = False)
    index_bed(f'{args.data_path}{out_name}')
    old_bed = fqtl_bed(normcount, f'{args.data_path}{args.gtf_file}')
    old_bed.to_csv(f'{args.data_path}old_{out_name}', sep = '\t', index = False)
    index_bed(f'{args.data_path}old_{out_name}')
    os.system(f'Rscript CRC_Generation.R {args.alias_out}.txt {args.cov_file} {args.ref_course}')

