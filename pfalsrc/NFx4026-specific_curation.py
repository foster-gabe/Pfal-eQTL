#!/usr/bin/env python
# coding: utf-8



######################################################################
#
# Name: NF54xNHP4026 Specific Data Curation
# Author: Gabe Foster
# Date: 6/1/2020
# Purpose: This script performs the initial preprocessing of the 
# NF54HT-GFP-luc x NHP4026 transcriptional cross section- this
# is the curation that is unique to this cross
#
######################################################################

import pandas as pd
import numpy as np
import re
import argparse
import os


def fix_names(countdata):
    '''
    This fixes the name specific issues in the NF54GFPxNHP4026 cross.
    The names have changed several times and been recorded in different
    formats, so I'll fix them with this.
    '''
    countdata = countdata.T
    countdata.index = countdata.index.str.replace('\/','', regex = True)
    countdata.index = countdata.index.str.replace('ND5A5', 'AC075', regex = True)
    countdata.index = countdata.index.str.replace('ND6G8', 'AC125', regex = True)
    countdata.index = countdata.index.str.replace('N1', '', regex = True)
    countdata.index = countdata.index.str.replace('\\.', '', regex = True)
    countdata.index = countdata.index.str.replace('_4026', '_NHP4026', regex = True)
    countdata.index = countdata.index.str.replace('^4026', 'NHP4026', regex = True)
    countdata.index = countdata.index.str.replace('2H9', 'AC030', regex = True)
    countdata.index = countdata.index.str.replace('6E5', 'AC033', regex = True)
    return countdata.T
   


def pick_best_reps(countdata):
    '''
    Get this, the GTEx Consortium handles TECHNICAL replicates
    (please don't mistake these for BIOLOGICAL replicates, which
    are amazing and useful) by simply choosing the replicate
    with the most reads. Good enough for the Broad, good enough
    for me, so I'll do that here. Note that this contains a
    pile of regular expressions for name handling specific to the
    formats used in this cross. 
    '''
    
    # Calculate sample count sums
    countdata = countdata.T
    good_samples = []
    sample_sum = []
    for i in range(6,len(countdata.index)):
        genes_rep = np.sum(countdata.iloc[i,:] > 0)
        good_samples.append(countdata.index[i])
        sample_sum.append(np.sum(countdata.iloc[i,:]))
        
    # Build summary frame, iterate over samples by strain/hpi and choose 
    # the sample with the most reads to keep

    count_summary = {'Sample Name':good_samples,
                     'Total Counts':sample_sum}
    count_summary = pd.DataFrame(count_summary)

    # This is a hot mess of regular expression that just converts the entire
    # sample name to strain_##, where ## is the sampling timepoint. The formatting
    # is very inconsistent throughout, so a mess of replacements need to be carefully
    # made. This mess of substitutions makes them.

    count_summary['Strain and Time'] = count_summary['Sample Name'].str.replace('^GF_PL[\d]+[a,b,c]{0,1}_', '',regex = True)
    count_summary['Strain and Time'] = count_summary['Strain and Time'].str.replace('[A,B,C]_', '', regex = True)
    count_summary['Strain and Time'] = count_summary['Strain and Time'].str.replace('_[d]+$', '', regex = True)
    count_summary['Strain and Time'] = count_summary['Strain and Time'].str.replace('_S.*', '', regex = True)
    count_summary['Strain and Time'] = count_summary['Strain and Time'].str.replace('_[0-9]{3,4}$', '', regex = True)
    count_summary['Strain and Time'] = count_summary['Strain and Time'].str.replace('hpi', '', regex = True)
    count_summary['Strain and Time'] = count_summary['Strain and Time'].str.replace('_T', '_', regex = True)

    # Believe it or not, the GTEx consortium handled reps across batches and such by...
    # taking the replicate with the most counts. That's easy, let's do that.

    best_samples = []
    for sample in count_summary['Strain and Time'].unique():
        subframe = count_summary[count_summary['Strain and Time'] == sample]
        best_samples.append(count_summary.iloc[subframe['Total Counts'].idxmax(),:])
    
    # Now I'll build a frame with just the best samples in it
    
    best_samples = pd.DataFrame(best_samples)
    best_samplenames = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
    best_samplenames.extend(list(best_samples['Sample Name']))
    curated_counts = countdata[countdata.index.isin(best_samplenames)]
    return curated_counts.T




def split_times(rawcounts):
    '''
    In this particular study I took 3 time points; 4hpi, 30hpi,
    and 44hpi. This function takes a full set of count data and
    splits it up in to 3 DataFrames. Note that it's hardcoded for the
    time points I took, so it's not suitable for all crosses.
    '''

    T4_samples = rawcounts.filter(regex='T4_|_4hpi')
    T4_samples = rawcounts.iloc[:,0:5].merge(T4_samples, 
                                             left_index = True, 
                                             right_index = True)
    T30_samples = rawcounts.filter(regex='T30_|_30hpi')
    T30_samples = rawcounts.iloc[:,0:5].merge(T30_samples, 
                                             left_index = True, 
                                             right_index = True)
    T44_samples = rawcounts.filter(regex='T44_|_44hpi')
    T44_samples = rawcounts.iloc[:,0:5].merge(T44_samples, 
                                             left_index = True, 
                                             right_index = True)



    # Let's pull out the count data for each time point
    # Now that we're normalized we can minimize this to a pure count matrix

    return [T4_samples, T30_samples, T44_samples]
    



def build_covariates(metadata):
    '''
    Our samples were run across different plates; this
    is the only real place I can account for technical
    variation. Our sampling batches are confounded with
    stage, so I am forced to rely on CRC and PEER to
    pull anything in there out. To run PEER, I need
    one-hot encoded plate data, so I do that here.
    '''
    encoding = pd.get_dummies(metadata['PlateID'])
    metadata = metadata.merge(encoding,
                              left_index = True, 
                              right_index = True)
    metadata.drop(columns = 'PlateID',
                 inplace = True)
    
    return metadata
    




def fix_vcf(vcf_file, progenydata, vcf_out):
    '''
    The vcf file contains TX versions of strain names; this
    function swaps them with the correct version.
    
    Additionally, it pulls only markers that differ between the parents. This is a simple
	thought that was a lot of effort.
    '''
    # Create dict for progeny data
    
    progenydata.index = progenydata.iloc[:,1]
    progenydict = progenydata.iloc[:,0]
    progenydict = progenydict.to_dict()

    # Swap names and write out- I tested using bcftools for this, and it was the same amount of time.

    vcf = open(vcf_file, 'r')
    vcfout = open(vcf_out, 'w')
    for line in vcf:
        for key in progenydict.keys():
            line = re.sub(key, progenydict[key], line)
        vcfout.write(line)
    
    vcf.close()
    vcfout.close()
	
	# Alright this is a mess- this block uses vcftools to create a matrix; this matrix is checked for
	# markers that differ between the parents, and only those markers are kept.

	# Pull handle for outputs
    vcfouthandle = vcf_out.replace('.vcf', '')
    
    # Write out genotype matrix
    os.system(f'vcftools --vcf {vcf_out} --012 --out {vcf_out}')
	
	# Read matrix, sample names and marker names - some name changing as well in here
    matrix = pd.read_csv(f'{vcf_out}.012', sep = '\t', header = None, index_col = 0)
    sample_names = pd.read_csv(f'{vcf_out}.012.indv', sep = '\t', header = None, index_col = False)
    sample_names.columns = ['samples']
    marker_names = pd.read_csv(f'{vcf_out}.012.pos', sep = '\t', header = None, index_col = False)
    marker_names.columns = ['CHROM', 'POS']
    marker_names['POS'] = marker_names['POS'].astype(str)
	
	# Generate marker IDs in the format CHROM~POS
    marker_names['ID'] = marker_names['CHROM'].str.cat(marker_names['POS'], sep = '~')
	
	#Write CHR POS ID frame out so we can annotate the new IDs to the VCF file
    marker_names.to_csv(f'{vcf_out}.anno', sep = '\t', header = False, index = False)
	
	# Annotate and index
    os.system(f'bgzip {vcf_out}.anno')
    os.system(f'tabix -p vcf {vcf_out}.anno.gz')
    os.system(f'bcftools annotate -c CHROM,POS,ID -a {vcf_out}.anno.gz {vcf_out} -o {vcfouthandle}_anno.vcf')
	
	# Here we just manipulate the matrix to just get the NF54 and 4026 rows and get the marker IDs
	# for markers that are different
	
    matrix.index = sample_names['samples']
    matrix.columns = marker_names['ID']

    dif_markers = matrix.loc['NHP4026',:] == matrix.loc['NF54gfp', :]
    dif_markers = dif_markers[dif_markers == False]
    dif_markers = dif_markers.index.to_series()
	
	# Write out list of markers to keep, and use vcftools to build a vcf with those only
	
    dif_markers.to_csv(f'{vcf_out}.diff', sep = '\t', header = False, index = False)
    os.system(f'vcftools --vcf {vcfouthandle}_anno.vcf --snps {vcf_out}.diff --recode --out {vcfouthandle}')
	
	# Clean up intermediate files
    os.system('rm /./data/*anno* /./data/*.log /./data/*.diff /./data/*012*')

    # finalmatrix = matrix.drop(dif_markers[dif_markers == True].index, axis = 1)
    # finalmatrix.to_csv(f'{vcf_out}.matrix', sep = '\t')
    

def strip_names(countdata):
    '''
    This function strips the long sample names from the
    count data down to their strain only. NOTE: You really
    don't want to use this until after you've built the metadata
    file, as you need the full sample names to get the correct
    Plate Number for batch correction.
    '''
    
    countdata.columns = countdata.columns.str.replace('^GF[\d]*_', '', regex = True)
    countdata.columns = countdata.columns.str.replace('PL[\da-z]*_', '', regex = True)
    countdata.columns = countdata.columns.str.replace('[A-C]{1}_', '', regex = True)
    countdata.columns = countdata.columns.str.partition('_').to_frame().iloc[:,0]
    
    return countdata

def to_gct(file, expdata):
	row = len(expdata.index)
	col = len(expdata.columns)
	outfile = open(file, 'w')
	outfile.write(f'1.0\n')
	outfile.write(f'{row}\t{col}\n')
	outfile.close()
	expdata.to_csv(file,
                       mode = 'a',
                       header = True,
                       sep = '\t',
                       encoding='utf-8',
                       index = False)

########################################################
#
# Main block- parses all files in one pass, if you're
# feeling lazy
#
########################################################




if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Curates data from the NF54GFP x NHP4026 Transcriptional Experiment for eQTL Analysis')
    parser.add_argument('counts_file', help = 'the *.count file provided to us after mapping')
    parser.add_argument('metadata', help = 'the *_name.txt file profixed to us after mapping')
    parser.add_argument('vcf_file', help = 'the *.vcf file from this genetic cross')
    parser.add_argument('vcf_out', help = 'file name for renamed vcf output')
    parser.add_argument('data_path', help = 'path for specific xls file for this cross only')
    args = parser.parse_args()

    
    # Read counts file, fix names, pick best reps, write out
    # probe metadata, and split out time points
    
    counts = pd.read_csv(f'{args.data_path}{args.counts_file}', sep = '\t')
    counts = fix_names(counts)
    counts = pick_best_reps(counts)
    
    probeinfo = counts.iloc[:,0:6]
    probeinfo.to_csv(f'{args.data_path}{args.counts_file.split(".")[0]}_probeinfo.csv')
    

 
    
    # Read in and correct metadata
    
    metadata = pd.read_csv(f'{args.data_path}{args.metadata}', sep = '\t', header = None, index_col = 0)
    metadata.columns = ['PlateID', 'Strain', 'Sampling Time', 'Sample Number', 'No Clue', 'Same']
    covariates = build_covariates(metadata)
    covariates.drop(columns = ['Strain', 'Sampling Time', 'Sample Number', 'No Clue', 'Same'], inplace = True)
    covariates = fix_names(covariates.T).T
	
    # build correct, separate metadata files for time points,
    # and write out metadata and expression
	
    timepoints = split_times(counts)
    
    timepoint_times = ['T4_', 'T30_', 'T44_']
    i = 0
    for timepoint in timepoints:
      metasub = covariates[covariates.index.isin(timepoint.columns)].T
      metasub = strip_names(metasub).T
	  # NOTE: we have RNAseq for this parasite, but no genome seq
      metasub.drop(index = 'AC081', inplace = True)
      metasub.to_csv(f'{args.data_path}{timepoint_times[i]}batchcov.csv')
      timepoint.drop(columns = ['Start', 'End', 'Strand', 'Chr'], inplace = True)
      timepoint = strip_names(timepoint)
      timepoint.drop(columns = 'AC081', inplace = True)
      timepoint.to_csv(f'{args.data_path}{timepoint_times[i]}counts.txt',
        sep = '\t',
        index = False)  
      i = i + 1
     
    # Fix names in vcf file
    
    progenydata = pd.read_excel(f'{args.data_path}NF54gfpluc NHP4026 progeny in map 6-25-18 .xlsx',
                            usecols = ['freezerPro ID', 'Map 7/5/18'])
    progenydata = progenydata.drop(progenydata.index[57])
    progenydata.iloc[55,0] = 'NHP4026'
    progenydata.iloc[56,0] = 'NF54gfp'
    fix_vcf(f'{args.data_path}{args.vcf_file}', progenydata, f'{args.data_path}{args.vcf_out}')
    




